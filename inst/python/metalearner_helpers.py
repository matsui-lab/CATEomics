"""
Helpers for extracting potential outcomes from fitted econml metalearners /
DR-Learner.

All helpers return (y0, y1) numpy arrays of length ``n = X.shape[0]``.

Design notes
------------
* T-Learner and X-Learner expose ``models`` (list of length 2). For these,
  ``models[0].predict(X)`` and ``models[1].predict(X)`` are the control /
  treated outcome surfaces, respectively.
* DR-Learner's ``models_regression`` is a nested list of shape
  ``[mc_iters][cv_folds]``. Each inner model takes ``[X | onehot(T)]`` as
  input. We compute μ_t(X) as the fold-wise average of predictions made with
  the T column set to ``t``. No re-fitting.
* This module only uses sklearn conventions and econml public attributes. It
  does not reach into private APIs beyond ``transformer`` / ``models_*``
  which are documented in EconML 0.15+.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np


def _as_2d(X: np.ndarray) -> np.ndarray:
    X = np.asarray(X)
    if X.ndim == 1:
        X = X.reshape(-1, 1)
    return X


def _mu_predict(model, X: np.ndarray) -> np.ndarray:
    """Predict E[Y|X,T=t] from a fold nuisance model.

    For a discrete (binary) outcome the fold model is a classifier, so return
    ``P(Y=1|X)`` from ``predict_proba`` (the column for class 1) rather than the
    hard 0/1 label ``predict`` would give; for a continuous outcome the
    regressor's ``predict`` is used directly.
    """
    if hasattr(model, "predict_proba"):
        proba = np.asarray(model.predict_proba(X))
        classes = getattr(model, "classes_", None)
        if classes is not None and 1 in list(classes):
            j = list(classes).index(1)
        else:
            j = proba.shape[1] - 1
        return proba[:, j].ravel()
    return np.asarray(model.predict(X)).ravel()


def _drlearner_augment(est, X: np.ndarray, t_value: int) -> np.ndarray:
    """Build DRLearner-style regression features ``[X | onehot(T)]``.

    DRLearner's ``models_regression`` are trained on ``[X | est.transformer(T)]``
    (no extra control column).
    """
    n = X.shape[0]
    t_col = np.full((n, 1), float(t_value))
    t_enc = est.transformer.transform(t_col)
    return np.concatenate([X, t_enc], axis=1)




def tlearner_potential_outcomes(est, X) -> Tuple[np.ndarray, np.ndarray]:
    """Return (y0, y1) for a fitted TLearner via est.models[0/1].predict(X)."""
    X = _as_2d(X)
    y0 = np.asarray(est.models[0].predict(X)).ravel()
    y1 = np.asarray(est.models[1].predict(X)).ravel()
    return y0, y1


def xlearner_potential_outcomes(est, X) -> Tuple[np.ndarray, np.ndarray]:
    """Return (y0, y1) for a fitted XLearner using Stage 1 est.models[0/1].

    mu0_hat, mu1_hat come from the Stage 1 outcome regressions (est.models),
    not the Stage 2 CATE models.
    """
    X = _as_2d(X)
    y0 = np.asarray(est.models[0].predict(X)).ravel()
    y1 = np.asarray(est.models[1].predict(X)).ravel()
    return y0, y1


def drlearner_potential_outcomes(est, X) -> Tuple[np.ndarray, np.ndarray]:
    """Return (y0, y1) for a fitted DRLearner via direct nuisance predictions.

      mu0_hat = mean over (mc_iters, folds) of models_regression.predict(X|T=0)
      mu1_hat = mean over (mc_iters, folds) of models_regression.predict(X|T=1)

    The orthogonalized CATE (``est.effect``) is used for ``ite`` separately; it
    is NOT required that ``mu1_hat - mu0_hat == ite`` for DR-Learner.
    """
    X = _as_2d(X)
    X0 = _drlearner_augment(est, X, 0)
    X1 = _drlearner_augment(est, X, 1)

    mr = est.models_regression  # list[mc_iters] of list[cv_folds]
    preds_0, preds_1 = [], []
    for mc_models in mr:
        for fold_model in mc_models:
            preds_0.append(_mu_predict(fold_model, X0))
            preds_1.append(_mu_predict(fold_model, X1))
    y0 = np.mean(np.vstack(preds_0), axis=0)
    y1 = np.mean(np.vstack(preds_1), axis=0)
    return y0, y1


def _propensity_predict(model, X: np.ndarray) -> np.ndarray:
    """Return P(T=1|X) from a fitted binary propensity classifier.

    econml reads the treated column as ``predict_proba(X)[:, 1:]``; look the
    class up by label so a reordered ``classes_`` cannot silently flip the
    score.
    """
    proba = np.asarray(model.predict_proba(X))
    classes = getattr(model, "classes_", None)
    if classes is not None and 1 in list(classes):
        j = list(classes).index(1)
    else:
        j = proba.shape[1] - 1
    return proba[:, j].ravel()


def xlearner_propensity(est, X) -> np.ndarray:
    """Return P(T=1|X) from a fitted XLearner.

    ``est.propensity_models`` holds one model per non-control treatment; with a
    binary treatment that is a single model fit on the whole sample, and it is
    the same object ``const_marginal_effect`` weights the two imputed-effect
    surfaces with.
    """
    X = _as_2d(X)
    models = getattr(est, "propensity_models", None)
    if not models:
        return None
    return _propensity_predict(models[0], X)


def drlearner_propensity(est, X) -> np.ndarray:
    """Return P(T=1|X) from a fitted DRLearner, averaged over the fold models.

    ``models_propensity`` is nested ``[mc_iters][cv_folds]``. econml does not
    persist which rows each fold held out (``fitted_inds`` is local to
    ``_crossfit``), so a training row cannot be matched to the one model that
    did not see it; the average over folds is what can be reconstructed. On
    rows the estimator never saw, every fold model is out-of-sample and the
    average is an ordinary ensemble prediction.
    """
    X = _as_2d(X)
    mp = getattr(est, "models_propensity", None)
    if not mp:
        return None
    preds = []
    for mc_models in mp:
        for fold_model in mc_models:
            preds.append(_propensity_predict(fold_model, X))
    if not preds:
        return None
    return np.mean(np.vstack(preds), axis=0)
