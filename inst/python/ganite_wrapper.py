"""cfomics wrapper for PyPI ganite 0.1.2 (Yoon 2018, BSD-3-Clause).

Wraps the PyPI `ganite` 0.1.2 (Yoon et al. official PyTorch implementation)
behind a cfomics-convention `run_ganite_pypi()`.

Backend characteristics:
- Backend: PyTorch
- GPU: works on Blackwell sm_120 natively (no PTX JIT 30+ min hang)
- z-score: PyPI ganite has NO internal normalization; X is z-scored centrally
  on the R side and Y is z-scored by methods_ganite.R
- RNG: 3-source save/restore (numpy + torch + random) around fit to prevent
  global state leak. PyPI ganite uses numpy global RNG for minibatch sampling
  (`np.random.randint`) and torch global RNG for weight initialization;
  `enable_reproducible_results()` helper is seed=0 hardcoded so cfomics sets
  the 3 RNGs explicitly per `random_state`.
"""

import random

import numpy as np
import torch
from ganite import Ganite


def run_ganite_pypi(
    X,
    T,
    Y,
    edges_list=None,           # kept for signature compat (GANITE ignores DAG)
    variable_names=None,       # kept for signature compat (provenance only)
    random_state=0,
    h_dim=100,
    alpha=2.0,
    beta=5.0,
    iterations=5000,
    batch_size=256,
    verbose=False,
    device="cpu",
):
    """Fit PyPI ganite 0.1.2 with isolated RNG and return cfomics result.

    Parameters
    ----------
    X : numpy.ndarray (n, p), float64
        Z-scored covariates (z-score performed on R side).
    T : numpy.ndarray (n,), float64
        Binary treatment (0/1).
    Y : numpy.ndarray (n,), float64
        Z-scored outcome.
    random_state : int
        Seed for numpy / torch / python-random global RNG (save/restore
        wraps the fit). Same seed = bit-identical fit.
    h_dim, alpha, beta, iterations, batch_size :
        Standard GANITE hyperparameters.
    device : str
        "cpu" or "gpu" or "cuda". "gpu" maps to "cuda" if available, else
        falls back to "cpu" with a stderr message. PyTorch backend supports
        Blackwell sm_120 natively.

    Returns
    -------
    dict
        {"model": Ganite instance, "res": {"mu0_hat", "mu1_hat"}}, where
        mu0_hat / mu1_hat are the PyPI ganite InferenceNets per-row outputs
        in the z-scored Y scale (caller denormalizes).
    """
    # ----- device resolution -----
    # NOTE: PyPI ganite (`ganite/ganite.py:12`) uses a module-global
    #   DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # and `.to(DEVICE)`s tensors inside its training loop unconditionally.
    # If the caller requests device="cpu" but cuda IS available, the inner
    # Ganite ops will still move things to cuda → mismatch with our cpu
    # input tensors. We override the module-global DEVICE here to respect
    # the caller's device choice.
    import ganite.ganite as _ganite_mod  # access to DEVICE constant
    if device in ("gpu", "cuda"):
        if torch.cuda.is_available():
            torch_device = torch.device("cuda")
        else:
            torch_device = torch.device("cpu")
            if verbose:
                print("ganite_wrapper: cuda requested but not available, "
                      "falling back to cpu.")
    else:
        torch_device = torch.device("cpu")
    # Save and override the upstream module-global DEVICE
    _saved_device = _ganite_mod.DEVICE
    _ganite_mod.DEVICE = torch_device

    # ----- RNG state save (3 sources: np + torch + random) -----
    np_state = np.random.get_state()
    torch_state = torch.get_rng_state()
    if torch.cuda.is_available():
        cuda_state = torch.cuda.get_rng_state_all()
    else:
        cuda_state = None
    random_state_py = random.getstate()

    try:
        # ----- seed all 3 RNGs -----
        # Ganite.__init__ calls enable_reproducible_results(), which upstream
        # hard-sets seed 0 and would overwrite any external seeding, making
        # random_state a no-op. Patch it (in this ephemeral worker process only)
        # to seed with random_state instead, so cf_fit's random_state actually
        # varies the ganite fit (weight init + minibatch sampling).
        _rs = int(random_state)

        def _seed_with_random_state():
            np.random.seed(_rs)
            torch.manual_seed(_rs)
            if torch.cuda.is_available():
                torch.cuda.manual_seed_all(_rs)
            random.seed(_rs)

        _ganite_mod.enable_reproducible_results = _seed_with_random_state
        _seed_with_random_state()

        # ----- prepare tensors on device -----
        X_arr = np.asarray(X, dtype=np.float32)
        T_arr = np.asarray(T, dtype=np.float32).reshape(-1, 1)
        Y_arr = np.asarray(Y, dtype=np.float32).reshape(-1, 1)
        X_t = torch.from_numpy(X_arr).to(torch_device)
        T_t = torch.from_numpy(T_arr).to(torch_device)
        Y_t = torch.from_numpy(Y_arr).to(torch_device)

        # ----- fit (Ganite.__init__ runs _fit internally) -----
        model = Ganite(
            X=X_t,
            Treatments=T_t,
            Y=Y_t,
            dim_hidden=int(h_dim),
            alpha=float(alpha),
            beta=float(beta),
            minibatch_size=int(batch_size),
            num_iterations=int(iterations),
        )

        # ----- predict mu0/mu1 separately via inference_nets -----
        # Ganite.forward returns ite = mu1 - mu0; we need both columns.
        # inference_nets(X) -> (n, 2) tensor with columns [mu0, mu1].
        with torch.no_grad():
            y_hat = model.inference_nets(X_t).detach().cpu().numpy()
        mu0_hat = y_hat[:, 0].astype(np.float64)
        mu1_hat = y_hat[:, 1].astype(np.float64)

        # Probe the actual device of inference_nets parameters so cf_fit can
        # record `device_effective` even when the DEVICE override silently
        # fell back.
        try:
            device_effective = str(next(model.inference_nets.parameters()).device)
        except StopIteration:
            device_effective = str(torch_device)

        return {
            "model": model,
            "res": {
                "mu0_hat": mu0_hat,
                "mu1_hat": mu1_hat,
                "device_effective": device_effective,
            },
        }
    finally:
        # ----- RNG state restore -----
        np.random.set_state(np_state)
        torch.set_rng_state(torch_state)
        if cuda_state is not None:
            torch.cuda.set_rng_state_all(cuda_state)
        random.setstate(random_state_py)
        # ----- DEVICE override restore -----
        _ganite_mod.DEVICE = _saved_device


def predict_ganite_pypi(model, X_new, device=None):
    """Predict mu0/mu1 on newdata via a fitted PyPI Ganite instance.

    methods_ganite.R `predict_cf_ganite` calls this for newdata.

    Auto-detects the device from `model.inference_nets` parameters when
    `device` is None. The caller
    (cf_fit / .kfold_oof_predict_matrix) may invoke predict on a model fit
    with `device="cuda"`; passing `device="cpu"` (the prior default) would
    create a cpu X_t while the model still holds cuda weights, raising
    "Expected all tensors to be on the same device" at the first
    `nn.Linear.forward` inside `inference_nets`. Auto-detection removes
    the caller's responsibility to track which device the fit used.
    """
    if device is None:
        # Probe one parameter from inference_nets to find the model's device.
        try:
            torch_device = next(model.inference_nets.parameters()).device
        except StopIteration:
            torch_device = torch.device("cpu")
    elif device in ("gpu", "cuda") and torch.cuda.is_available():
        torch_device = torch.device("cuda")
    elif device in ("gpu", "cuda"):
        torch_device = torch.device("cpu")
    else:
        torch_device = torch.device(device)

    X_arr = np.asarray(X_new, dtype=np.float32)
    X_t = torch.from_numpy(X_arr).to(torch_device)
    with torch.no_grad():
        y_hat = model.inference_nets(X_t).detach().cpu().numpy()
    mu0_hat = y_hat[:, 0].astype(np.float64)
    mu1_hat = y_hat[:, 1].astype(np.float64)
    return {"mu0_hat": mu0_hat, "mu1_hat": mu1_hat}


def save_ganite(model, path):
    """Serialise a fitted GANITE inference network to `path`.

    Only ``model.inference_nets`` is needed to score new samples, so the
    adversarial training sub-networks are not persisted.
    """
    import torch
    torch.save(model.inference_nets, path)


class _ReloadedGanite:
    """Minimal object exposing ``.inference_nets`` for predict_ganite_pypi."""

    def __init__(self, inference_nets):
        self.inference_nets = inference_nets


def load_ganite(path, device=None):
    """Reload a GANITE inference network saved by :func:`save_ganite`."""
    import torch
    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    infnets = torch.load(path, map_location=device, weights_only=False)
    return _ReloadedGanite(infnets)
