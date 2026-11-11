"""CATEomics wrapper for catenets TARNet / DragonNet.

This thin module isolates the reticulate <-> Python boundary so the
CATEomics R side can:

  * fit a TARNet / DragonNet via :func:`fit_tarnet` / :func:`fit_dragonnet`
  * receive `(cate, mu0_hat, mu1_hat)` as ``np.ndarray`` triplets
  * deterministically control the RNG with a single ``random_state`` integer

`catenets v0.2.x` keeps the ``seed`` constructor argument for
train/validation splitting, but the per-batch shuffling and weight
initialisation use the *global* numpy / torch / Python RNG state. This
wrapper saves and restores those globals around each fit, so the seeding
stays local to the fit.

TARNet (Shalit-Johansson-Sontag ICML 2017) and DragonNet (Shi-Blei-Veitch
NeurIPS 2019) are provided via catenets v0.2.x. CFRNet semantics are obtained
from TARNet with ``penalty_disc > 0``.
"""

from __future__ import annotations

import contextlib
import random as _py_random
import sys
import numpy as np
import torch

from catenets.models.torch import TARNet, DragonNet


@contextlib.contextmanager
def _pinned_device(device):
    """Run the body with catenets' torch device forced to ``device``.

    catenets v0.2.x takes no device argument. Each of its torch modules
    holds a module-level ``DEVICE = torch.device("cuda" if
    torch.cuda.is_available() else "cpu")``, decided when the module is
    first imported, and its training loops call ``.to(DEVICE)``
    unconditionally. Hiding the GPU with ``CUDA_VISIBLE_DEVICES`` from the
    R side cannot steer that: the Python work runs in a basilisk worker
    process that is reused across calls, and a process inherits the
    environment only at spawn time, so a variable set afterwards never
    reaches an interpreter that is already running. Rebind the constant
    instead, and restore it after. catenets is the only backend that needs
    this: the others take the device as an argument on every call.

    ``sys.modules`` is scanned rather than a fixed list of module names
    because several catenets modules take their own copy of the constant
    with ``from ... import DEVICE``, so rebinding the defining module alone
    would leave those copies pointing at the GPU.
    """
    want = str(device).lower()
    if want in ("gpu", "cuda", "auto"):
        target = (torch.device("cuda") if torch.cuda.is_available()
                  else torch.device("cpu"))
    elif want in ("cpu", "", "none"):
        target = torch.device("cpu")
    else:
        raise ValueError("device must be one of cpu, gpu, cuda, auto; got %r"
                         % (device,))
    holders = [
        m for name, m in list(sys.modules.items())
        if name.startswith("catenets") and m is not None
        and getattr(m, "DEVICE", None) is not None
    ]
    saved = [(m, m.DEVICE) for m in holders]
    for m in holders:
        m.DEVICE = target
    try:
        yield target
    finally:
        for m, old in saved:
            m.DEVICE = old


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
    val_split_prop=0.3,
    device="cpu",
    **kwargs,
):
    """Fit a catenets TARNet.

    Returns a dict with keys ``model``, ``mu0_hat``, ``mu1_hat``, ``ite``.

    catenets v0.2.x does not expose a direct `device` parameter, so the
    requested device is applied by :func:`_pinned_device`, which rebinds the
    module-level ``DEVICE`` constant catenets reads. ``device='gpu'`` falls
    back to the CPU when no GPU is visible.
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
            # The fraction catenets holds out to decide when to stop; the
            # weights see the rest. Explicit here so the default is visible
            # rather than inherited silently from catenets.
            val_split_prop=float(val_split_prop),
        )
        if random_state is not None:
            ctor_kwargs["seed"] = int(random_state)
        # Extra hyperparameters come from the caller's `params` block, so an
        # unknown name is a user error: let TARNet() raise TypeError.
        ctor_kwargs.update(kwargs)

        with _pinned_device(device):
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
    val_split_prop=0.3,
    device="cpu",  # see fit_tarnet docstring for device handling
    **kwargs,
):
    """Fit a catenets DragonNet.

    DragonNet (Shi-Blei-Veitch 2019) shares the TARNet representation
    head and adds an auxiliary propensity head, whose binary cross-entropy
    enters the training loss weighted by ``prop_loss_multiplier``. catenets
    v0.2.x implements no targeted regularisation.

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
            # The fraction catenets holds out to decide when to stop; the
            # weights see the rest. Explicit here so the default is visible
            # rather than inherited silently from catenets.
            val_split_prop=float(val_split_prop),
        )
        if random_state is not None:
            ctor_kwargs["seed"] = int(random_state)
        # See fit_tarnet: an unknown name is a user error, so DragonNet()
        # raises rather than the wrapper dropping it.
        ctor_kwargs.update(kwargs)

        with _pinned_device(device):
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


def predict_with(model, X, device="cpu"):
    """Predict (cate, mu0, mu1) on new X for a previously-fit model.

    Pinned the same way as the fit: the forward pass also reads the
    module-level ``DEVICE`` to place its input tensors, so leaving it alone
    would move a CPU-trained model's prediction onto whatever GPU happens to
    be visible.
    """
    X_np = np.asarray(X, dtype=np.float32)
    with _pinned_device(device):
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
    if device is None or (str(device).startswith("cuda")
                          and not torch.cuda.is_available()):
        # A model fitted on a GPU can be asked for on a host without one.
        device = "cuda" if torch.cuda.is_available() else "cpu"
    model = torch.load(path, map_location=device, weights_only=False)
    try:
        model.device = torch.device(device)
    except Exception:
        pass
    return model
