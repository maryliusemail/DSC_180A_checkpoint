"""Covariate-adjusted eQTL regression and genotype-PCA helpers."""
import numpy as np
from scipy import stats


def covariate_basis(covariates):
    C = np.asarray(covariates, dtype=float)
    if not np.isfinite(C).all() or np.linalg.matrix_rank(C) != C.shape[1]:
        raise ValueError("Covariates must be finite and have independent columns.")
    return np.linalg.qr(C, mode="reduced")[0]


def associate_adjusted(y, X, indices, covariates, cache=None):
    """Two-sided partial OLS tests, with exact pairwise missingness and residual df."""
    indices = np.asarray(indices, dtype=int)
    if not len(indices):
        return tuple(np.array([], dtype=int if i < 2 else float) for i in range(5))
    if cache is not None and np.isfinite(y).all():
        Q = cache["Q"]
        yr = y - Q @ (Q.T @ y)
        yy = yr @ yr
        if yy <= 1e-12:
            return tuple(np.array([], dtype=int if i < 2 else float) for i in range(5))
        complete = cache["complete"][indices]
        good_indices = indices[complete & (cache["xx"][indices] > 1e-12)]
        pieces = []
        if len(good_indices):
            xx = cache["xx"][good_indices]
            xy = cache["R"][:, good_indices].T @ yr
            beta = xy / xx
            df = len(y) - covariates.shape[1] - 1
            se = np.sqrt(np.maximum(yy - xy * beta, 0) / df / xx)
            with np.errstate(divide="ignore", invalid="ignore"):
                p = 2 * stats.t.sf(np.abs(beta / se), df)
            valid = np.isfinite(beta) & np.isfinite(se) & np.isfinite(p)
            pieces.append((good_indices[valid], np.full(valid.sum(), len(y), dtype=int), beta[valid], se[valid], p[valid]))
        if (~complete).any():
            pieces.append(associate_adjusted(y, X, indices[~complete], covariates))
        if not pieces:
            return tuple(np.array([], dtype=int if i < 2 else float) for i in range(5))
        values = tuple(np.concatenate([p[i] for p in pieces]) for i in range(5))
        order = np.argsort(values[0])
        return tuple(v[order] for v in values)
    x = X[:, indices]
    ok = np.isfinite(x) & np.isfinite(y[:, None]) & np.isfinite(covariates).all(axis=1)[:, None]
    if ok.all():
        groups = [(np.ones(len(y), dtype=bool), np.arange(len(indices)))]
    else:
        _, labels = np.unique(np.packbits(ok, axis=0).T, axis=0, return_inverse=True)
        groups = [(ok[:, np.flatnonzero(labels == g)[0]], np.flatnonzero(labels == g))
                  for g in np.unique(labels)]
    pieces = []
    for rows, columns in groups:
        C = covariates[rows]
        n = int(rows.sum())
        if n <= C.shape[1] + 1 or np.linalg.matrix_rank(C) != C.shape[1]:
            continue
        Q = covariate_basis(C)
        yr = y[rows] - Q @ (Q.T @ y[rows])
        yy = yr @ yr
        if yy <= 1e-12:
            continue
        xr = x[np.ix_(rows, columns)]
        xr = xr - Q @ (Q.T @ xr)
        xx = np.einsum("ij,ij->j", xr, xr)
        valid = xx > 1e-12
        if not valid.any():
            continue
        ix = indices[columns[valid]]
        xy = xr[:, valid].T @ yr
        beta = xy / xx[valid]
        df = n - C.shape[1] - 1
        se = np.sqrt(np.maximum(yy - xy * beta, 0) / df / xx[valid])
        with np.errstate(divide="ignore", invalid="ignore"):
            p = 2 * stats.t.sf(np.abs(beta / se), df)
        finite = np.isfinite(beta) & np.isfinite(se) & np.isfinite(p)
        pieces.append((ix[finite], np.full(finite.sum(), n, dtype=int), beta[finite], se[finite], p[finite]))
    if not pieces:
        return tuple(np.array([], dtype=int if i < 2 else float) for i in range(5))
    values = tuple(np.concatenate([p[i] for p in pieces]) for i in range(5))
    order = np.argsort(values[0])
    return tuple(v[order] for v in values)


def check_adjusted_pair(y, dosage, covariates, beta, se, p):
    complete = np.isfinite(y) & np.isfinite(dosage) & np.isfinite(covariates).all(axis=1)
    D = np.column_stack([covariates[complete], dosage[complete]])
    yy = y[complete]
    coefficients = np.linalg.lstsq(D, yy, rcond=None)[0]
    df = len(yy) - D.shape[1]
    direct_se = np.sqrt(np.sum((yy - D @ coefficients)**2) / df * np.linalg.inv(D.T @ D)[-1, -1])
    direct_p = 2 * stats.t.sf(abs(coefficients[-1] / direct_se), df)
    np.testing.assert_allclose([beta, se, p], [coefficients[-1], direct_se, direct_p], rtol=1e-6, atol=1e-10)


def relationship_context(X, positions, marker_mask):
    """Standardized, pruned genotype markers for gene-window-excluded PCA."""
    Z = X[:, marker_mask].copy()
    f = Z.mean(axis=0) / 2
    if not np.isfinite(Z).all() or not ((f > 0) & (f < 1)).all():
        raise ValueError("PCA markers must be complete and variable.")
    Z -= 2 * f
    Z /= np.sqrt(2 * f * (1 - f))
    return {"Z": Z, "positions": np.asarray(positions)[marker_mask], "K": Z @ Z.T}


def population_covariates(context, start, end, is_yri, n_pcs=5):
    """Joint chromosome 22 PCs, excluding the current gene's entire test window."""
    inside = (context["positions"] >= start) & (context["positions"] <= end)
    remaining = context["Z"].shape[1] - int(inside.sum())
    if remaining < 100:
        raise ValueError("Too few PCA markers outside the cis window.")
    removed = context["Z"][:, inside]
    K = (context["K"] - removed @ removed.T) / remaining
    values, vectors = np.linalg.eigh(K)
    pcs = vectors[:, np.argsort(values)[::-1][:n_pcs]].copy()
    # Population intercepts are already in the model. Centering within groups
    # preserves the model's column space and improves numerical conditioning.
    for group in (False, True):
        rows = np.asarray(is_yri, dtype=bool) == group
        pcs[rows] -= pcs[rows].mean(axis=0)
    scale = pcs.std(axis=0)
    if (scale < 1e-10).any():
        raise ValueError("A genotype PC is redundant with population membership.")
    pcs /= scale
    C = np.column_stack([np.ones(len(is_yri)), np.asarray(is_yri, float), pcs])
    covariate_basis(C)
    return C, remaining


def interaction_hc3(y, X, is_yri, covariates):
    """SNP + SNP:population, conditional on intercept/population/PC covariates.

    Returns SNP effect in EUR, effect in YRI, difference, HC3 SE, p-value,
    and mask identifying estimable SNPs. Wald inference uses the normal limit.
    """
    Q = covariate_basis(covariates)
    yr = y - Q @ (Q.T @ y)
    U = X - Q @ (Q.T @ X)
    interaction = X * np.asarray(is_yri, float)[:, None]
    V = interaction - Q @ (Q.T @ interaction)
    a = np.einsum("ij,ij->j", U, U)
    b = np.einsum("ij,ij->j", U, V)
    c = np.einsum("ij,ij->j", V, V)
    determinant = a * c - b * b
    estimable = (a > 1e-10) & (c > 1e-10) & (determinant > 1e-10 * np.maximum(a * c, 1))
    safe = np.where(estimable, determinant, 1)
    uy, vy = U.T @ yr, V.T @ yr
    beta_eur = (c * uy - b * vy) / safe
    difference = (a * vy - b * uy) / safe
    residual = yr[:, None] - U * beta_eur - V * difference
    h0 = np.einsum("ij,ij->i", Q, Q)[:, None]
    leverage = h0 + (c * U * U - 2 * b * U * V + a * V * V) / safe
    influence = (-b * U + a * V) / safe
    with np.errstate(divide="ignore", invalid="ignore"):
        variance = ((influence * residual / (1 - leverage))**2).sum(axis=0)
        se = np.sqrt(variance)
        p = 2 * stats.norm.sf(np.abs(difference / se))
    valid = estimable & np.isfinite(p) & np.isfinite(se) & (se > 0) & (leverage < 1 - 1e-10).all(axis=0)
    return beta_eur[valid], (beta_eur + difference)[valid], difference[valid], se[valid], p[valid], valid


def check_interaction(y, dosage, is_yri, covariates, difference, se, p):
    D = np.column_stack([covariates, dosage, dosage * np.asarray(is_yri, float)])
    inverse = np.linalg.inv(D.T @ D)
    coefficients = np.linalg.lstsq(D, y, rcond=None)[0]
    residual = y - D @ coefficients
    leverage = np.einsum("ij,jk,ik->i", D, inverse, D)
    weighted = D * (residual / (1 - leverage))[:, None]
    covariance = inverse @ (weighted.T @ weighted) @ inverse
    direct_se = np.sqrt(covariance[-1, -1])
    direct_p = 2 * stats.norm.sf(abs(coefficients[-1] / direct_se))
    np.testing.assert_allclose([difference, se, p], [coefficients[-1], direct_se, direct_p], rtol=1e-5, atol=1e-9)


def genotype_adjustment_cache(X, covariates):
    """Residualize complete-call SNPs once per chromosome, not once per gene."""
    Q = covariate_basis(covariates)
    complete = np.isfinite(X).all(axis=0)
    R = np.zeros_like(X, dtype=float)
    R[:, complete] = X[:, complete] - Q @ (Q.T @ X[:, complete])
    return {"Q": Q, "complete": complete, "R": R, "xx": np.einsum("ij,ij->j", R, R)}
