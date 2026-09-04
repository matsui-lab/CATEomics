"""cfomics wrapper for catenets TARNet / DragonNet.

This thin module isolates the reticulate <-> Python boundary so the
cfomics R side can:

  * fit a TARNet / DragonNet via :func:`fit_tarnet` / :func:`fit_dragonnet`
  * receive `(cate, mu0_hat, mu1_hat)` as ``np.ndarray`` triplets
  * deterministically control the RNG with a single ``random_state`` integer

`catenets v0.2.x` keeps the ``seed`` constructor argument for
train/validation splitting, but the per-batch shuffling and weight
initialisation use the *global* numpy / torch / Python RNG state. This
wrapper saves and restores those globals around each fit so cfomics tests
get bit-exact reproducibility (max diff = 0).

TARNet (Shalit-Johansson-Sontag ICML 2017) and DragonNet (Shi-Blei-Veitch
NeurIPS 2019) are provided via catenets v0.2.x. CFRNet semantics are obtained
from TARNet with ``penalty_disc > 0``.
"""

from __future__ import annotations

import random as _py_random
import numpy as np
import torch

from catenets.models.torch import TARNet, DragonNet


def _set_global_seeds(seed: int) -> None:
    """Force the three RNGs catenets touches at fit time."""
    s = int(seed)
    np.random.seed(s)
    torch.manual_seed(s)
    _py_random.seed(s)


def _save_rng_state() -> dict:
    return {
        "np": np.random.get_state(),
        "torch": torch.get_rng_state(),
        "py": _py_random.getstate(),
    }


def _restore_rng_state(state: dict) -> None:
    np.random.set_state(state["np"])
    torch.set_rng_state(state["torch"])
    _py_random.setstate(state["py"])


def _to_np(t):
    """Convert a torch tensor (or anything tensor-like) to a 1-D np.ndarray."""
    if hasattr(t, "detach"):
        t = t.detach().cpu().numpy()
    arr = np.asarray(t).reshape(-1)
    return arr.astype(np.float64, copy=False)


def _get_effective_device(module):
    """Return the str device of the first parameter of a torch nn.Module.

    Callers receive the *actual* device the model ended up on (e.g.
    ``"cuda:0"`` or ``"cpu"``), rather than the device that was *requested*
    by the user. This lets cf_fit record both `device_requested` and
    `device_effective` in `meta$provenance` so downstream provenance audits
    are accurate even when an upstream library silently falls back (e.g.
    CUDA libs missing).
    """
    try:
        return str(next(module.parameters()).device)
    except StopIteration:
        return "cpu"


def fit_tarnet(
    X,
    T,
    Y,
    random_state=None,
    binary_y=False,
    penalty_disc=0.0,
    n_iter=10000,
    batch_size=100,
    learning_rate=1e-4,
    device="cpu",
    **kwargs,
):
    """Fit a catenets TARNet.

    Returns a dict with keys ``model``, ``mu0_hat``, ``mu1_hat``, ``ite``.

    catenets v0.2.x does not expose a direct `device` parameter, so GPU
    enable/disable is controlled on the cfomics side via the
    `CUDA_VISIBLE_DEVICES` env var (cf_fit.R `.resolve_device()` sets it
    empty when `device='cpu'`). With `device='gpu'` and a GPU available,
    catenets' internal PyTorch uses CUDA automatically.
    """
    X_np = np.asarray(X, dtype=np.float32)
    T_np = np.asarray(T, dtype=np.float32)
    Y_np = np.asarray(Y, dtype=np.float32)

    saved = _save_rng_state()
    try:
        if random_state is not None:
            _set_global_seeds(random_state)

        ctor_kwargs = dict(
            n_unit_in=X_np.shape[1],
            binary_y=bool(binary_y),
            penalty_disc=float(penalty_disc),
            n_iter=int(n_iter),
            batch_size=int(batch_size),
            lr=float(learning_rate),
        )
        if random_state is not None:
            ctor_kwargs["seed"] = int(random_state)
        # Filter kwargs by TARNet constructor signature. cf_compare forwards
        # `...` to every method's fit, so per-method args from siblings (e.g.
        # ganite `iterations`, cevae `num_epochs`) reach here via **kwargs and
        # TARNet() rejects them with TypeError. Same pattern as the R-side
        # `formals()` filter in cf_fit_grf.
        import inspect as _inspect
        _allowed = set(_inspect.signature(TARNet.__init__).parameters)
        ctor_kwargs.update({k: v for k, v in kwargs.items() if k in _allowed})

        model = TARNet(**ctor_kwargs)
        model.fit(X_np, Y_np, T_np)

        cate_t, mu0_t, mu1_t = model.predict(X_np, return_po=True)
        cate = _to_np(cate_t)
        mu0 = _to_np(mu0_t)
        mu1 = _to_np(mu1_t)
        device_effective = _get_effective_device(model)
    finally:
        _restore_rng_state(saved)

    return {
        "model": model,
        "mu0_hat": mu0,
        "mu1_hat": mu1,
        "ite": cate,
        "device_effective": device_effective,
    }


def fit_dragonnet(
    X,
    T,
    Y,
    random_state=None,
    binary_y=False,
    n_iter=10000,
    batch_size=100,
    learning_rate=1e-4,
    device="cpu",  # see fit_tarnet docstring for device handling
    **kwargs,
):
    """Fit a catenets DragonNet.

    DragonNet (Shi-Blei-Veitch 2019) shares the TARNet representation
    head with an auxiliary propensity head; in catenets v0.2.x the
    targeted-regularisation option is enabled internally when
    ``binary_y=False``.

    Returns a dict with the same shape as :func:`fit_tarnet`.
    """
    X_np = np.asarray(X, dtype=np.float32)
    T_np = np.asarray(T, dtype=np.float32)
    Y_np = np.asarray(Y, dtype=np.float32)

    saved = _save_rng_state()
    try:
        if random_state is not None:
            _set_global_seeds(random_state)

        ctor_kwargs = dict(
            n_unit_in=X_np.shape[1],
            binary_y=bool(binary_y),
            n_iter=int(n_iter),
            batch_size=int(batch_size),
            lr=float(learning_rate),
        )
        if random_state is not None:
            ctor_kwargs["seed"] = int(random_state)
        # Filter kwargs by DragonNet ctor signature; see fit_tarnet for the
        # same handling.
        import inspect as _inspect
        _allowed = set(_inspect.signature(DragonNet.__init__).parameters)
        ctor_kwargs.update({k: v for k, v in kwargs.items() if k in _allowed})

        model = DragonNet(**ctor_kwargs)
        model.fit(X_np, Y_np, T_np)

        cate_t, mu0_t, mu1_t = model.predict(X_np, return_po=True)
        cate = _to_np(cate_t)
        mu0 = _to_np(mu0_t)
        mu1 = _to_np(mu1_t)
        device_effective = _get_effective_device(model)
    finally:
        _restore_rng_state(saved)

    return {
        "model": model,
        "mu0_hat": mu0,
        "mu1_hat": mu1,
        "ite": cate,
        "device_effective": device_effective,
    }


def predict_with(model, X):
    """Predict (cate, mu0, mu1) on new X for a previously-fit model."""
    X_np = np.asarray(X, dtype=np.float32)
    cate_t, mu0_t, mu1_t = model.predict(X_np, return_po=True)
    return {
        "ite": _to_np(cate_t),
        "mu0_hat": _to_np(mu0_t),
        "mu1_hat": _to_np(mu1_t),
    }


def save_model(model, path):
    """Serialise a fitted catenets model to `path`."""
    import torch
    torch.save(model, path)


def load_model(path, device=None):
    """Reload a catenets model saved by :func:`save_model`.

    Loads onto CUDA when available, otherwise CPU. A whole-catenets model
    saved from one device needs its parameter tensors *and* its internal
    ``.device`` attribute to agree, or ``predict`` raises a device mismatch,
    so both are set to the resolved device.
    """
    import torch
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    model = torch.load(path, map_location=device, weights_only=False)
    try:
        model.device = torch.device(device)
    except Exception:
        pass
    return model
