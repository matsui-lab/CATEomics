"""GANITE (Yoon, Jordon and van der Schaar, ICLR 2018) in PyTorch.

Four sub-networks, in two blocks:

  counterfactual block   G    generator, trunk + two outcome heads
                         D_G  tells which of the two outcomes was observed
  ITE block              I    inference network, trunk + two outcome heads
                         D_I  tells G's completed outcomes from I's

The blocks train one after the other: G and D_G first, then I and D_I on the
dataset G completed. Each block runs `iterations` steps, and in each step the
discriminator is updated `d_steps` times before its counterpart once.

The heads emit a value on the outcome scale. With a binary outcome they emit a
logit and every consumer -- the discriminator input, the targets the ITE block
regresses, the predictions -- reads the probability, so the supervised factual
term is a cross-entropy while the ITE term stays a squared difference of
probabilities (a risk difference has no cross-entropy).
"""

import numpy as np
import torch
import torch.nn as nn


def _hidden(in_dim, width, depth):
    """`depth` hidden ReLU layers of `width`, and nothing after them."""
    layers = []
    d = in_dim
    for _ in range(max(1, int(depth))):
        layers += [nn.Linear(d, width), nn.ReLU()]
        d = width
    return nn.Sequential(*layers)


def _mlp(in_dim, width, depth, out_dim):
    """`depth` hidden ReLU layers of `width`, then a linear output."""
    return nn.Sequential(_hidden(in_dim, width, depth),
                         nn.Linear(width, out_dim))


def _xavier(module):
    """Xavier uniform on every linear layer; torch defaults to Kaiming."""
    for m in module.modules():
        if isinstance(m, nn.Linear):
            nn.init.xavier_uniform_(m.weight)
            nn.init.zeros_(m.bias)
    return module


class _TwoHead(nn.Module):
    """Trunk with one head per arm, as in the paper's multi-task layout."""

    def __init__(self, in_dim, width, depth):
        super().__init__()
        # `depth` layers in the trunk and one more in each head.
        self.trunk = _hidden(in_dim, width, depth)
        self.head0 = _mlp(width, width, 1, 1)
        self.head1 = _mlp(width, width, 1, 1)

    def forward(self, x):
        h = self.trunk(x)
        return self.head0(h), self.head1(h)


class GaniteModel:
    """The four networks, their optimisers and the two-block training loop."""

    def __init__(self, input_dim, h_dim=100, num_layers=2, alpha=2.0,
                 beta=5.0, learning_rate=1e-3, weight_decay=1e-4,
                 d_steps=2, binary_y=False, device="cpu", random_state=0):
        self.input_dim = int(input_dim)
        self.h_dim = int(h_dim)
        self.num_layers = int(num_layers)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.d_steps = max(1, int(d_steps))
        self.binary_y = bool(binary_y)
        self.device = torch.device(device)
        self.random_state = int(random_state)

        torch.manual_seed(self.random_state)
        if self.device.type == "cuda":
            torch.cuda.manual_seed_all(self.random_state)

        p = self.input_dim
        # G sees (x, t, y_factual); D_G and D_I see (x, y0, y1).
        self.G = _xavier(_TwoHead(p + 2, self.h_dim, self.num_layers)).to(self.device)
        self.I = _xavier(_TwoHead(p, self.h_dim, self.num_layers)).to(self.device)
        self.D_G = _xavier(_mlp(p + 2, self.h_dim, self.num_layers, 1)).to(self.device)
        self.D_I = _xavier(_mlp(p + 2, self.h_dim, self.num_layers, 1)).to(self.device)

        opt = lambda m: torch.optim.Adam(m.parameters(), lr=float(learning_rate),
                                         weight_decay=float(weight_decay))
        self.opt_G, self.opt_I = opt(self.G), opt(self.I)
        self.opt_DG, self.opt_DI = opt(self.D_G), opt(self.D_I)
        self.bce = nn.BCEWithLogitsLoss()
        self.training_trace = []

    # -- scales ----------------------------------------------------------
    def _value(self, head_out):
        """What every consumer reads: a probability when the outcome is binary."""
        return torch.sigmoid(head_out) if self.binary_y else head_out

    # -- forwards --------------------------------------------------------
    def _g_heads(self, x, t, y):
        return self.G(torch.cat([x, t, y], dim=1))

    def _complete(self, x, t, y, detach=False):
        """G's completed pair: the observed outcome on the assigned arm."""
        g0, g1 = self._g_heads(x, t, y)
        v0, v1 = self._value(g0), self._value(g1)
        if detach:
            v0, v1 = v0.detach(), v1.detach()
        y0 = (1.0 - t) * y + t * v0
        y1 = t * y + (1.0 - t) * v1
        return y0, y1, g0, g1

    def _dg_logit(self, x, y0, y1):
        return self.D_G(torch.cat([x, y0, y1], dim=1))

    def _di_logit(self, x, y0, y1):
        return self.D_I(torch.cat([x, y0, y1], dim=1))

    def _i_values(self, x):
        i0, i1 = self.I(x)
        return self._value(i0), self._value(i1)

    # -- steps -----------------------------------------------------------
    def _step_DG(self, x, t, y):
        """D_G: which of the pair was observed, i.e. the treatment label."""
        y0, y1, _, _ = self._complete(x, t, y, detach=True)
        loss = self.bce(self._dg_logit(x, y0, y1), t)
        self.opt_DG.zero_grad(set_to_none=True)
        loss.backward()
        self.opt_DG.step()
        return float(loss.detach())

    def _step_G(self, x, t, y):
        """G: fool D_G, and match the outcome that was observed."""
        y0, y1, g0, g1 = self._complete(x, t, y)
        # Saturating form, which is what `min_G V_CF` says: G minimises the
        # negative of what D_G minimises. The ITE block uses the non-saturating
        # form instead because there the constant term is the one that carries
        # D_I's own pair, so nothing would reach I (see `_step_I`).
        adv = -self.bce(self._dg_logit(x, y0, y1), t)
        if self.binary_y:
            logit_factual = t * g1 + (1.0 - t) * g0
            factual = self.bce(logit_factual, y)
        else:
            value_factual = t * self._value(g1) + (1.0 - t) * self._value(g0)
            factual = torch.mean((y - value_factual) ** 2)
        loss = adv + self.alpha * factual
        self.opt_G.zero_grad(set_to_none=True)
        loss.backward()
        self.opt_G.step()
        return float(loss.detach()), float(factual.detach())

    def _step_DI(self, x, t, y):
        """D_I: G's completed pair (real) against I's pair (fake).

        The returned loss is the sum of the two cross-entropies, so chance
        level is 2 log 2 rather than the log 2 of `_step_DG`: the `d_g` and
        `d_i` columns of the trace are not on one scale.
        """
        with torch.no_grad():
            y0_bar, y1_bar, _, _ = self._complete(x, t, y, detach=True)
            i0, i1 = self._i_values(x)
        real = self._di_logit(x, y0_bar, y1_bar)
        fake = self._di_logit(x, i0, i1)
        loss = (self.bce(real, torch.ones_like(real)) +
                self.bce(fake, torch.zeros_like(fake)))
        self.opt_DI.zero_grad(set_to_none=True)
        loss.backward()
        self.opt_DI.step()
        return float(loss.detach())

    def _step_I(self, x, t, y):
        """I: have D_I take its pair for G's, and match G's effect."""
        with torch.no_grad():
            y0_bar, y1_bar, _, _ = self._complete(x, t, y, detach=True)
        i0, i1 = self._i_values(x)
        fake = self._di_logit(x, i0, i1)
        # Non-saturating form (Goodfellow et al. 2014): the term that reaches I
        # is the one that asks D_I to call its pair real.
        adv = self.bce(fake, torch.ones_like(fake))
        l_s_i = torch.mean(((y1_bar - y0_bar) - (i1 - i0)) ** 2)
        loss = adv + self.beta * l_s_i
        self.opt_I.zero_grad(set_to_none=True)
        loss.backward()
        self.opt_I.step()
        return float(loss.detach()), float(l_s_i.detach())

    # -- fit / predict ---------------------------------------------------
    def fit(self, X, T, Y, iterations=1000, batch_size=256, trace_every=100):
        rng = np.random.default_rng(self.random_state)
        x_all = torch.as_tensor(np.asarray(X, dtype=np.float32), device=self.device)
        t_all = torch.as_tensor(np.asarray(T, dtype=np.float32).reshape(-1, 1),
                                device=self.device)
        y_all = torch.as_tensor(np.asarray(Y, dtype=np.float32).reshape(-1, 1),
                                device=self.device)
        n = x_all.shape[0]
        size = min(int(batch_size), n)

        def batch():
            idx = rng.choice(n, size=size, replace=False)
            i = torch.as_tensor(idx, device=self.device)
            return x_all[i], t_all[i], y_all[i]

        iterations = max(1, int(iterations))
        # `iterations` counts the steps of each block, so the two blocks run
        # `iterations` each. Every `trace_every`-th step is recorded, and the
        # last one always is, so the trace shows where training ended.
        for it in range(iterations):
            for _ in range(self.d_steps):
                d_g = self._step_DG(*batch())
            g, g_factual = self._step_G(*batch())
            if it % trace_every == 0 or it == iterations - 1:
                self.training_trace.append(
                    dict(block="counterfactual", iteration=it, d_g=d_g, g=g,
                         g_factual=g_factual, d_i=float("nan"), i=float("nan"),
                         l_s_i=float("nan")))
        for it in range(iterations):
            for _ in range(self.d_steps):
                d_i = self._step_DI(*batch())
            i, l_s_i = self._step_I(*batch())
            if it % trace_every == 0 or it == iterations - 1:
                self.training_trace.append(
                    dict(block="ite", iteration=it, d_g=float("nan"),
                         g=float("nan"), g_factual=float("nan"), d_i=d_i, i=i,
                         l_s_i=l_s_i))
        return self

    @torch.no_grad()
    def predict(self, X):
        """(mu0, mu1) from the ITE block, so mu1 - mu0 is the effect exactly."""
        x = torch.as_tensor(np.asarray(X, dtype=np.float32), device=self.device)
        i0, i1 = self._i_values(x)
        return (i0.cpu().numpy().ravel(), i1.cpu().numpy().ravel())

    # -- serialisation ---------------------------------------------------
    def ite_state(self):
        """Everything predict() needs: the ITE block and how to rebuild it."""
        return dict(
            config=dict(input_dim=self.input_dim, h_dim=self.h_dim,
                        num_layers=self.num_layers, binary_y=self.binary_y),
            state_dict={k: v.cpu() for k, v in self.I.state_dict().items()})


class GanitePredictor:
    """The ITE block on its own, rebuilt from a saved state."""

    def __init__(self, state, device="cpu"):
        cfg = state["config"]
        self.binary_y = bool(cfg["binary_y"])
        self.device = torch.device(device)
        self.I = _TwoHead(int(cfg["input_dim"]), int(cfg["h_dim"]),
                          int(cfg["num_layers"])).to(self.device)
        self.I.load_state_dict(state["state_dict"])
        self.I.eval()

    @torch.no_grad()
    def predict(self, X):
        x = torch.as_tensor(np.asarray(X, dtype=np.float32), device=self.device)
        i0, i1 = self.I(x)
        if self.binary_y:
            i0, i1 = torch.sigmoid(i0), torch.sigmoid(i1)
        return (i0.cpu().numpy().ravel(), i1.cpu().numpy().ravel())
