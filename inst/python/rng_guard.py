"""Save and restore the Python-side global RNG state around a CATEomics call.

basilisk evaluates the worker body in the calling R process whenever its
shared mode is active, so ``numpy``, ``torch`` and ``random`` keep their global
state between calls. Every fitter seeds those generators from ``random_state``,
which would otherwise leave the caller's own Python session advanced by a
CATEomics fit. The R-side generator is already saved and restored by the methods
that call ``set.seed()``; these two helpers do the same for the Python side.

``snapshot()`` is cheap (three small objects, plus one per visible CUDA
device) and tolerates a missing torch, so it can wrap every call without a
per-method branch.
"""

import random as _py_random

import numpy as np


def snapshot():
    """Return the current global RNG state of numpy / random / torch."""
    state = {
        "np": np.random.get_state(),
        "py": _py_random.getstate(),
    }
    try:
        import torch
    except Exception:
        return state
    state["torch"] = torch.get_rng_state()
    if torch.cuda.is_available():
        state["cuda"] = torch.cuda.get_rng_state_all()
    return state


def restore(state):
    """Put back a state returned by :func:`snapshot`."""
    if not state:
        return
    if "np" in state:
        np.random.set_state(state["np"])
    if "py" in state:
        _py_random.setstate(state["py"])
    if "torch" in state:
        import torch

        torch.set_rng_state(state["torch"])
        if "cuda" in state and torch.cuda.is_available():
            torch.cuda.set_rng_state_all(state["cuda"])
