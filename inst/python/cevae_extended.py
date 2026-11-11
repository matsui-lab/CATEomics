"""CEVAE extension exposing individual y0 and y1 predictions.

Pyro's built-in `CEVAE.ite()` returns only the difference `y1 - y0`, hiding the
two potential outcomes. CATEomics needs both individually for the `mu0_hat` /
`mu1_hat` contract. `CEVAEExtended.predict_y(x)` mirrors the internals of the
upstream `ite()` (py-why/pyro `pyro/contrib/cevae/__init__.py`, commit as of
pyro 1.9.1) but returns `(y0_mean, y1_mean)` instead of collapsing them.

`CEVAEExtended.fit()` overrides the upstream `fit()` method because the
upstream constructs `torch.Generator(device=x.device)` with an OS-sourced
seed for DataLoader shuffling. That generator ignores `torch.manual_seed()`,
so same-seed repeat runs otherwise produce different minibatch orders and
therefore different parameter updates. We inject a generator seeded from the
torch global RNG (which `cf_fit_cevae` seeds deterministically) so repeated
calls with the same `random_state` are reproducible.
"""

import random

import numpy as np
import pyro
import pyro.distributions as dist
import pyro.poutine as poutine
import torch
from pyro.contrib.cevae import CEVAE, TraceCausalEffect_ELBO
from pyro.infer import SVI
from pyro.optim import ClippedAdam
from torch.utils.data import DataLoader, TensorDataset


class _DeviceAwareZDist:
    """Picklable device-aware latent-z prior for CEVAE.

    A module-level class (not a ``fit``-local closure) so the whole model can
    be serialised with ``torch.save`` (a local function is not picklable). The
    ``device`` attribute is re-synced by :func:`load_cevae` after reload so
    the prior matches wherever the parameters were mapped.
    """

    def __init__(self, device, latent_dim):
        self.device = device
        self.latent_dim = latent_dim

    def __call__(self):
        loc = torch.zeros(self.latent_dim, device=self.device)
        scale = torch.ones(self.latent_dim, device=self.device)
        return dist.Normal(loc, scale).to_event(1)


class _NoWhiten(torch.nn.Module):
    """Identity replacement for pyro's ``PreWhitener``.

    CATEomics standardises the covariate matrix centrally on the R side
    (``.cf_preprocess_input``: every column z-scored, binary indicators
    included). Pyro's ``PreWhitener`` would standardise again on top of that,
    so it is replaced by this no-op.
    Kept as a picklable ``nn.Module`` so the fitted model still serialises with
    ``torch.save``.
    """

    def forward(self, x):
        return x


class CEVAEExtended(CEVAE):
    """CEVAE subclass that exposes per-sample `y0` and `y1` expectations."""

    def get_effective_device(self):
        """Return str of the actual device the model's parameters live on.

        cf_fit reads this to populate `meta$provenance$device_effective`, so
        audits see the device the model really ended up on (not just what was
        requested).
        """
        try:
            return str(next(self.parameters()).device)
        except StopIteration:
            return "cpu"

    def fit(
        self,
        x,
        t,
        y,
        num_epochs=100,
        batch_size=100,
        learning_rate=1e-3,
        learning_rate_decay=0.1,
        weight_decay=1e-4,
        log_every=100,
    ):
        """Seed-aware override of `pyro.contrib.cevae.CEVAE.fit`.

        Follows the upstream implementation (pyro 1.9.1) with these
        differences: the `torch.Generator` used for DataLoader shuffling is
        seeded from the torch global RNG (upstream passes
        `torch.Generator(device=x.device)` without a seed, so the shuffle
        order is not reproducible; the seed is derived via `torch.randint`
        against the already-seeded global generator); the parameter store is
        cleared and the module moved to the data's device before fitting;
        `PreWhitener` is replaced with `_NoWhiten`, because the design matrix
        arrives standardised and whitening it twice changes the numbers; the
        latent prior is taken from the override of `model.z_dist`; and
        upstream's assertion that the loss is not NaN is not repeated, so a
        run that diverges returns rather than stopping inside pyro.
        """
        assert x.dim() == 2 and x.size(-1) == self.feature_dim
        assert t.shape == x.shape[:1]
        assert y.shape == y.shape[:1]

        # Pyro's `pyro.module("cevae", ...)` registers params in the global
        # store, so a second `fit` on a freshly-constructed CEVAEExtended
        # would otherwise inherit the first fit's trained parameters under
        # the same name. Clear the store so each fit starts from the
        # (seeded) newly-initialised module parameters.
        pyro.clear_param_store()

        # Move the entire CEVAE Module subtree to x.device so that all internal
        # nn parameters (Normal.loc, BernoulliNet params, etc.) align with the
        # data tensors. Without this, the constructor's cpu-default initial
        # weights mismatch a cuda x_tensor during svi.step → "Expected all
        # tensors to be on the same device".
        self.to(x.device)

        # X is standardised centrally on the R side, so the pyro PreWhitener is
        # replaced by an identity (see `_NoWhiten`) to avoid standardising
        # twice. `.to(x.device)` is a no-op for the buffer-less
        # identity but kept for symmetry.
        self.whiten = _NoWhiten().to(x.device)

        # Override `Model.z_dist` to make the latent z prior device-aware.
        # Upstream pyro 1.9.1 returns
        # `dist.Normal(0, 1).expand([latent_dim]).to_event(1)`; the scalar
        # 0/1 create cpu tensors regardless of where `self.to(device)` moved
        # the nn parameters. The model.forward then samples z from a cpu
        # Normal while x_dist/y_dist net outputs live on cuda, and the SVI
        # log_prob computation hits `Expected all tensors to be on the same
        # device, cuda:0 and cpu`. The override binds the z prior to the
        # caller's device. Set on the instance once and inherited by both
        # fit (svi.step) and predict (y_mean) calls thereafter.
        _device = x.device
        _latent_dim = self.model.latent_dim

        # Picklable callable (module-level class) so the fitted model can be
        # serialised with torch.save; a fit-local closure could not.
        self.model.z_dist = _DeviceAwareZDist(_device, _latent_dim)
        # Cache device on model so predict_y can rebuild device-aware do() tensors
        self.model._cateomics_device = _device

        # PyTorch 2.x DataLoader internally seeds via a cpu int64 tensor:
        #   self._base_seed = torch.empty((), dtype=torch.int64).random_(
        #       generator=loader.generator).item()
        # That destination tensor is cpu regardless of where the data tensors
        # live, so `loader.generator` must be cpu. Upstream pyro CEVAE uses
        # `torch.Generator(device=x.device)`, which fires "Expected a 'cpu'
        # device type for generator but found 'cuda'" when x is on cuda.
        # Force a cpu generator here — the seed source only feeds DataLoader's
        # cpu-side shuffle bookkeeping, so it does not need to match x.device.
        loader_seed = int(torch.randint(0, 2**31 - 1, (1,)).item())
        loader_generator = torch.Generator(device="cpu")
        loader_generator.manual_seed(loader_seed)

        dataset = TensorDataset(x, t, y)
        dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=True,
            generator=loader_generator,
        )
        num_steps = num_epochs * len(dataloader)
        optim = ClippedAdam(
            {
                "lr": learning_rate,
                "weight_decay": weight_decay,
                "lrd": learning_rate_decay ** (1 / num_steps),
            }
        )
        svi = SVI(self.model, self.guide, optim, TraceCausalEffect_ELBO())
        losses = []
        for epoch in range(num_epochs):
            for xb, tb, yb in dataloader:
                xb = self.whiten(xb)
                loss = svi.step(xb, tb, yb, size=len(dataset)) / len(dataset)
                losses.append(loss)
        return losses

    @torch.no_grad()
    def predict_y(self, x, num_samples=None, batch_size=None):
        """Return the Monte-Carlo averaged y0 and y1 for each row of `x`.

        Mirrors `CEVAE.ite()` control flow (including optional batching and
        the feature-dim assertion) so that `predict_y(x)[1] - predict_y(x)[0]`
        equals the upstream `ite(x)` up to Monte-Carlo noise from a separate
        sample of z.
        """
        if num_samples is None:
            num_samples = self.num_samples
        if not torch._C._get_tracing_state():
            assert x.dim() == 2 and x.size(-1) == self.feature_dim

        # Rebuild t={0,1} interventions on x.device. Upstream pyro CEVAE uses
        # bare `torch.zeros(())` / `torch.ones(())` (cpu 0-dim scalars) which
        # mismatch when x is on cuda (the model's z_dist has been overridden to
        # x.device in fit, and the model nn params are on cuda). Match the
        # scalar device to x so y_dist(t=., z=cuda) ops do not raise device
        # mismatch.
        _dev = x.device
        _t_zero = torch.zeros((), device=_dev)
        _t_one = torch.ones((), device=_dev)

        dataloader = [x] if batch_size is None else DataLoader(x, batch_size=batch_size)
        y0_chunks = []
        y1_chunks = []
        for xb in dataloader:
            xb = self.whiten(xb)
            with pyro.plate("num_particles", num_samples, dim=-2):
                with poutine.trace() as tr, poutine.block(hide=["y", "t"]):
                    self.guide(xb)
                with poutine.do(data=dict(t=_t_zero)):
                    y0 = poutine.replay(self.model.y_mean, tr.trace)(xb)
                with poutine.do(data=dict(t=_t_one)):
                    y1 = poutine.replay(self.model.y_mean, tr.trace)(xb)
            y0_chunks.append(y0.mean(0))
            y1_chunks.append(y1.mean(0))
        return torch.cat(y0_chunks), torch.cat(y1_chunks)

    @torch.no_grad()
    def predict_y_seeded(self, x, seed, num_samples=None, batch_size=None):
        """`predict_y` with a reproducible, isolated RNG.

        Saves the current torch / python-random / numpy RNG state, reseeds
        them from `seed`, runs `predict_y`, then restores the original state.
        This keeps repeat `predict(..., newdata=...)` calls deterministic for
        a given (fit, newdata, seed) tuple without leaking our seed into the
        caller's global RNG (which `pyro.set_rng_seed` would otherwise do).
        """
        t_state = torch.get_rng_state()
        p_state = random.getstate()
        n_state = np.random.get_state()
        try:
            torch.manual_seed(int(seed))
            random.seed(int(seed))
            np.random.seed(int(seed))
            return self.predict_y(
                x, num_samples=num_samples, batch_size=batch_size
            )
        finally:
            torch.set_rng_state(t_state)
            random.setstate(p_state)
            np.random.set_state(n_state)

    @torch.no_grad()
    def predict_y_with_ci(self, x, alpha=0.05, num_samples=None,
                          batch_size=None):
        """Return Pyro Predictive-style point estimates + CI for y0 / y1 / ite / ate.

        Instead of averaging the Monte-Carlo posterior samples before
        returning (as `predict_y` does), this variant retains the per-sample
        draws and computes percentile-method CIs over `num_samples` posterior
        samples of z. Output is a dict with eight fields (`y0_hat`, `y1_hat`,
        `ite`, `ate`, `ate_ci_lower`, `ate_ci_upper`, `ite_ci_lower`,
        `ite_ci_upper`) so the R side can populate `res$ate_ci` / `res$tau_ci`.
        """
        if num_samples is None:
            num_samples = self.num_samples
        if not torch._C._get_tracing_state():
            assert x.dim() == 2 and x.size(-1) == self.feature_dim

        # Same device-aware t={0,1} as predict_y.
        _dev = x.device
        _t_zero = torch.zeros((), device=_dev)
        _t_one = torch.ones((), device=_dev)

        dataloader = (
            [x] if batch_size is None else DataLoader(x, batch_size=batch_size)
        )
        y0_chunks = []
        y1_chunks = []
        for xb in dataloader:
            xb = self.whiten(xb)
            with pyro.plate("num_particles", num_samples, dim=-2):
                with poutine.trace() as tr, poutine.block(hide=["y", "t"]):
                    self.guide(xb)
                with poutine.do(data=dict(t=_t_zero)):
                    y0 = poutine.replay(self.model.y_mean, tr.trace)(xb)
                with poutine.do(data=dict(t=_t_one)):
                    y1 = poutine.replay(self.model.y_mean, tr.trace)(xb)
            # Keep the leading num_samples dimension for CI computation.
            y0_chunks.append(y0)
            y1_chunks.append(y1)
        # Each chunk shape: (num_samples, batch). Concat along the batch dim.
        y0_post = torch.cat(y0_chunks, dim=-1).cpu().numpy()
        y1_post = torch.cat(y1_chunks, dim=-1).cpu().numpy()
        ite_post = y1_post - y0_post

        y0_hat = y0_post.mean(axis=0)
        y1_hat = y1_post.mean(axis=0)
        ite = ite_post.mean(axis=0)
        ate = float(ite.mean())
        ate_post = ite_post.mean(axis=-1)
        lo = alpha / 2
        hi = 1 - alpha / 2
        return {
            "y0_hat": y0_hat,
            "y1_hat": y1_hat,
            "ite": ite,
            "ate": ate,
            "ate_ci_lower": float(np.quantile(ate_post, lo)),
            "ate_ci_upper": float(np.quantile(ate_post, hi)),
            "ite_ci_lower": np.quantile(ite_post, lo, axis=0),
            "ite_ci_upper": np.quantile(ite_post, hi, axis=0),
        }

    @torch.no_grad()
    def predict_y_with_ci_seeded(self, x, seed, alpha=0.05,
                                  num_samples=None, batch_size=None):
        """`predict_y_with_ci` with a reproducible, isolated RNG.

        Mirrors `predict_y_seeded` (saves / restores torch / random / numpy
        global RNG state around the call) so repeat `predict(..., newdata)`
        invocations on the same fit yield identical CIs without polluting
        the caller's global RNG.
        """
        t_state = torch.get_rng_state()
        p_state = random.getstate()
        n_state = np.random.get_state()
        try:
            torch.manual_seed(int(seed))
            random.seed(int(seed))
            np.random.seed(int(seed))
            return self.predict_y_with_ci(
                x, alpha=alpha, num_samples=num_samples, batch_size=batch_size
            )
        finally:
            torch.set_rng_state(t_state)
            random.setstate(p_state)
            np.random.set_state(n_state)

    def fit_seeded(self, x, t, y, seed, **kwargs):
        """`CEVAE.fit` with a reproducible, isolated RNG.

        Saves torch / python-random / numpy RNG state, reseeds them from
        `seed`, runs the upstream `fit`, then restores the original state.
        This makes `cf_fit(method="cevae", random_state=s)` deterministic
        for a given (x, t, y, s) tuple and prevents our seed from leaking
        into the caller's global RNG (which `pyro.set_rng_seed` — invoked
        internally by `CEVAE.fit` — would otherwise do).
        """
        t_state = torch.get_rng_state()
        p_state = random.getstate()
        n_state = np.random.get_state()
        try:
            torch.manual_seed(int(seed))
            random.seed(int(seed))
            np.random.seed(int(seed))
            return self.fit(x, t, y, **kwargs)
        finally:
            torch.set_rng_state(t_state)
            random.setstate(p_state)
            np.random.set_state(n_state)




def save_cevae(model, path):
    """Serialise a fitted CEVAEExtended to `path`."""
    import torch
    torch.save(model, path)


def load_cevae(path, device=None):
    """Reload a CEVAEExtended saved by :func:`save_cevae`.

    The nn-module parameters live on the model itself (predict reads the
    device from a registered parameter), so Pyro's global param store is not
    needed at predict time. The device-aware latent-z prior and the cached
    ``_cateomics_device`` are re-synced to wherever the parameters were mapped
    so the prior samples land on the same device as the parameters.
    """
    import torch
    if device in ("gpu", "cuda") or device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    model = torch.load(path, map_location=device, weights_only=False)
    try:
        dev = next(model.parameters()).device
        if hasattr(model.model, "z_dist") and hasattr(model.model.z_dist, "device"):
            model.model.z_dist.device = dev
        model.model._cateomics_device = dev
    except Exception:
        pass
    return model
