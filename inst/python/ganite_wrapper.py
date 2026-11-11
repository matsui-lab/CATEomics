"""CATEomics entry points for the GANITE implementation in `ganite_impl.py`.

Backend characteristics:
- Backend: PyTorch (`ganite_impl.GaniteModel`).
- z-score: the implementation normalises nothing; X is z-scored centrally on
  the R side and Y by `methods_ganite.R`.
- RNG: seeded from `random_state` inside the model (torch and a numpy
  Generator for the mini-batches). The caller's global numpy / random / torch
  state is snapshotted and restored by `.cf_py_run()`, so nothing leaks out.
- Serialisation: only the ITE block is saved, since that is what scores new
  samples.
"""

import numpy as np
import torch

from ganite_impl import GaniteModel, GanitePredictor


def _resolve_device(device):
    if device in ("gpu", "cuda"):
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device is None:
        return torch.device("cpu")
    return torch.device(device)


def run_ganite(X, T, Y, random_state=0, binary_y=False, h_dim=100,
               num_layers=2, alpha=2.0, beta=5.0, iterations=1000,
               batch_size=256, learning_rate=1e-3, weight_decay=1e-4,
               d_steps=2, device="cpu", **kwargs):
    """Fit GANITE and return the model plus its in-sample estimates.

    `iterations` counts the steps of each block, so the counterfactual block
    and the ITE block each run that many.

    Returns a dict with `model`, `mu0_hat`, `mu1_hat`, `ite`,
    `device_effective` and `training_trace`.
    """
    if kwargs:
        raise TypeError(
            "unexpected GANITE argument(s): " + ", ".join(sorted(kwargs)))

    X_arr = np.asarray(X, dtype=np.float32)
    if X_arr.ndim == 1:
        X_arr = X_arr.reshape(-1, 1)
    T_arr = np.asarray(T, dtype=np.float32).ravel()
    Y_arr = np.asarray(Y, dtype=np.float32).ravel()
    n = X_arr.shape[0]
    if T_arr.shape[0] != n or Y_arr.shape[0] != n:
        raise ValueError(
            f"shape mismatch: X has {n} rows, T has {T_arr.shape[0]}, "
            f"Y has {Y_arr.shape[0]}")
    uniq = np.unique(T_arr)
    if not np.all(np.isin(uniq, [0.0, 1.0])):
        raise ValueError(f"T must be binary (0/1), got: {uniq}")

    torch_device = _resolve_device(device)
    model = GaniteModel(
        input_dim=X_arr.shape[1], h_dim=h_dim, num_layers=num_layers,
        alpha=alpha, beta=beta, learning_rate=learning_rate,
        weight_decay=weight_decay, d_steps=d_steps, binary_y=binary_y,
        device=str(torch_device), random_state=random_state)
    model.fit(X_arr, T_arr, Y_arr, iterations=iterations,
              batch_size=batch_size)

    mu0, mu1 = model.predict(X_arr)
    try:
        device_effective = str(next(model.I.parameters()).device)
    except StopIteration:
        device_effective = str(torch_device)

    return {
        "model": model,
        "mu0_hat": mu0.astype(np.float64),
        "mu1_hat": mu1.astype(np.float64),
        "ite": (mu1 - mu0).astype(np.float64),
        "device_effective": device_effective,
        "training_trace": model.training_trace,
    }


def predict_ganite(model, X_new):
    """Score new samples with a fitted model or a reloaded ITE block.

    The model already carries its device, set when it was built or reloaded,
    so there is nothing to choose here.
    """
    X_arr = np.asarray(X_new, dtype=np.float32)
    if X_arr.ndim == 1:
        X_arr = X_arr.reshape(-1, 1)
    mu0, mu1 = model.predict(X_arr)
    return {"mu0_hat": mu0.astype(np.float64),
            "mu1_hat": mu1.astype(np.float64),
            "ite": (mu1 - mu0).astype(np.float64)}


def save_ganite(model, path):
    """Persist the ITE block; the adversarial networks are not needed again."""
    torch.save(model.ite_state(), path)


def load_ganite(path, device=None):
    """Rebuild the ITE block saved by :func:`save_ganite`."""
    torch_device = _resolve_device(device)
    state = torch.load(path, map_location=str(torch_device), weights_only=False)
    return GanitePredictor(state, device=str(torch_device))
