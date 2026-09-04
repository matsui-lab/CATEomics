"""sklearn-compatible wrapper that clips ``predict_proba`` to ``[L, U]``.

Supports the unified ``trim = c(L, U)`` API across PS-using methods.
``econml.dr.DRLearner`` exposes ``min_propensity`` (lower bound only) and
econml estimators expose ``trimming_threshold`` (symmetric only), so neither
natively supports the user-facing ``trim = c(L, U)`` asymmetric two-sided API.
We therefore wrap any sklearn-style probabilistic classifier (e.g., logistic /
lasso / RandomForestClassifier from cfomics' learner factory) and clip both
treated and control columns of ``predict_proba`` to ``[lower, upper]`` while
renormalising to keep rows summing to 1.

This file is imported lazily from R via ``reticulate::import_from_path`` so it
adds no Python load cost when ``trim`` is not in use.
"""

import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin, clone


class ClippedPropensityModel(ClassifierMixin, BaseEstimator):
    """Wrap a sklearn classifier and clip ``predict_proba`` output.

    Parameters
    ----------
    base_estimator : sklearn classifier
        Must implement ``fit(X, y)``, ``predict_proba(X)``, ``predict(X)``,
        and expose ``classes_`` after fitting.
    lower : float, default ``0.025``
        Lower clip threshold (Crump 2009 / van der Laan-Rose 2011 standard).
    upper : float, default ``0.975``
        Upper clip threshold.
    """

    # The (ClassifierMixin, BaseEstimator) MRO order (mixin first) is what lets
    # sklearn's tag system mark this wrapper as a classifier, so DoubleML does
    # not emit a spurious "not a classifier" warning. The explicit
    # `_estimator_type` is belt-and-suspenders for older sklearn.
    _estimator_type = "classifier"

    def __init__(self, base_estimator, lower=0.025, upper=0.975):
        if not (0.0 < lower < upper < 1.0):
            raise ValueError(
                f"Require 0 < lower < upper < 1; got lower={lower}, upper={upper}."
            )
        self.base_estimator = base_estimator
        self.lower = lower
        self.upper = upper

    def fit(self, X, y, sample_weight=None):
        self.estimator_ = clone(self.base_estimator)
        try:
            self.estimator_.fit(X, y, sample_weight=sample_weight)
        except TypeError:
            self.estimator_.fit(X, y)
        # Forward attributes for sklearn compatibility.
        self.classes_ = self.estimator_.classes_
        return self

    def predict_proba(self, X):
        proba = self.estimator_.predict_proba(X)
        proba_clipped = np.clip(proba, self.lower, self.upper)
        # Renormalise rows so they sum to 1.
        row_sums = proba_clipped.sum(axis=1, keepdims=True)
        # Guard against pathological all-zero rows (shouldn't happen after clip
        # given lower > 0, but defensively keep them as uniform if so).
        zero_mask = row_sums.flatten() == 0
        if np.any(zero_mask):
            proba_clipped[zero_mask] = 1.0 / proba_clipped.shape[1]
            row_sums = proba_clipped.sum(axis=1, keepdims=True)
        return proba_clipped / row_sums

    def predict(self, X):
        return self.estimator_.predict(X)

    def __repr__(self):
        return (
            f"ClippedPropensityModel(base_estimator={self.estimator_!r}, "
            f"lower={self.lower}, upper={self.upper})"
            if hasattr(self, "estimator_") else
            f"ClippedPropensityModel(base_estimator={self.base_estimator!r}, "
            f"lower={self.lower}, upper={self.upper})"
        )


class FixedPropensityModel(ClassifierMixin, BaseEstimator):
    """Return a known propensity (scalar or per-sample) regardless of X.

    Used by the known-PS path (RCT context) for econml methods that
    internally call ``predict_proba(X)``. ``predict_proba(X)`` returns a
    matrix where column 1 contains the requested ``propensity_score`` and
    column 0 is its complement, broadcasting a scalar to all rows.

    Parameters
    ----------
    propensity_score : float or 1d array
        If scalar, broadcast to all rows. If 1d array, length must match the
        number of rows seen in ``fit``.
    classes : tuple, default ``(0, 1)``
        Class labels reported via ``classes_``.
    """

    _estimator_type = "classifier"

    def __init__(self, propensity_score, classes=(0, 1)):
        self.propensity_score = propensity_score
        self.classes = classes

    def fit(self, X, y, sample_weight=None):
        # No fitting required; record the training size for vector validation.
        self.n_fit_ = X.shape[0]
        self.classes_ = np.array(self.classes)
        # Robustly detect scalar vs vector: reticulate often boxes R scalars
        # as length-1 numpy arrays, so size==1 also counts as scalar.
        ps = np.asarray(self.propensity_score, dtype=float).flatten()
        if ps.size == 1:
            self._ps_array_ = None
            self._scalar_ps_ = float(ps[0])
        elif ps.size == self.n_fit_:
            self._ps_array_ = ps
            self._scalar_ps_ = None
        else:
            raise ValueError(
                f"propensity_score must be a scalar or a 1d array of length n "
                f"(got size={ps.size}, n_fit={self.n_fit_})."
            )
        return self

    def predict_proba(self, X):
        n_new = X.shape[0]
        if self._ps_array_ is None:
            ps_col1 = np.full(n_new, self._scalar_ps_)
        else:
            if n_new != self.n_fit_:
                raise ValueError(
                    "predict_proba with newdata is unsupported for "
                    "vector-valued known propensity (length mismatch)."
                )
            ps_col1 = self._ps_array_
        proba = np.column_stack([1.0 - ps_col1, ps_col1])
        return proba

    def predict(self, X):
        proba = self.predict_proba(X)
        return self.classes_[(proba[:, 1] >= 0.5).astype(int)]


class ColumnSubsetPropensityModel(ClassifierMixin, BaseEstimator):
    """Wrap a probabilistic classifier so it only sees a column subset of X.

    The user can pass `propensity_X = c("X1", "X3")` to cf_fit / cf_compare to
    restrict PS estimation to a subset of covariates while the outcome model
    still uses the full X (DML standard practice; collider / mediator / proxy
    hygiene on the T-side). For sklearn-based methods (aipw / drlearner /
    xlearner) the
    natural way to express "fit the classifier on a column subset" is via
    `ColumnTransformer` + `Pipeline`, but that requires a multi-step glue per
    backend. This single-class wrapper achieves the same thing with no
    Pipeline plumbing: it slices X by `column_indices` before delegating
    `fit` / `predict` / `predict_proba` to the underlying estimator, then
    forwards `classes_` for sklearn-compat.

    Optionally combines column subsetting with the same `[lower, upper]`
    `predict_proba` clipping that `ClippedPropensityModel` provides — the
    PS-using method can therefore replace one wrapper with this combined
    wrapper (no double-wrapping).

    Parameters
    ----------
    base_estimator : sklearn classifier
    column_indices : 1d array-like of int
        0-based indices of the columns of X to keep before fitting. The
        unused columns are dropped from both `fit(X, y)` and
        `predict_proba(X)`.
    lower, upper : float, default ``0.025``, ``0.975``
        Optional clip thresholds (set both to ``None`` to disable clipping).
    """

    _estimator_type = "classifier"

    def __init__(self, base_estimator, column_indices,
                 lower=0.025, upper=0.975):
        if lower is not None and upper is not None:
            if not (0.0 < lower < upper < 1.0):
                raise ValueError(
                    f"Require 0 < lower < upper < 1; got lower={lower}, upper={upper}."
                )
        self.base_estimator = base_estimator
        self.column_indices = column_indices
        self.lower = lower
        self.upper = upper

    def _subset(self, X):
        idx = np.asarray(self.column_indices, dtype=int).flatten()
        if idx.ndim != 1:
            raise ValueError("column_indices must be a 1d array of int.")
        X_arr = np.asarray(X)
        if X_arr.ndim != 2:
            raise ValueError(f"X must be 2d; got shape {X_arr.shape}.")
        return X_arr[:, idx]

    def fit(self, X, y, sample_weight=None):
        X_sub = self._subset(X)
        self.estimator_ = clone(self.base_estimator)
        try:
            self.estimator_.fit(X_sub, y, sample_weight=sample_weight)
        except TypeError:
            self.estimator_.fit(X_sub, y)
        self.classes_ = self.estimator_.classes_
        return self

    def predict_proba(self, X):
        X_sub = self._subset(X)
        proba = self.estimator_.predict_proba(X_sub)
        if self.lower is None or self.upper is None:
            return proba
        proba_clipped = np.clip(proba, self.lower, self.upper)
        row_sums = proba_clipped.sum(axis=1, keepdims=True)
        zero_mask = row_sums.flatten() == 0
        if np.any(zero_mask):
            proba_clipped[zero_mask] = 1.0 / proba_clipped.shape[1]
            row_sums = proba_clipped.sum(axis=1, keepdims=True)
        return proba_clipped / row_sums

    def predict(self, X):
        return self.estimator_.predict(self._subset(X))

    def __repr__(self):
        n_cols = (len(self.column_indices)
                  if self.column_indices is not None else "?")
        return (
            f"ColumnSubsetPropensityModel(base_estimator={self.base_estimator!r}, "
            f"n_cols={n_cols}, lower={self.lower}, upper={self.upper})"
        )
