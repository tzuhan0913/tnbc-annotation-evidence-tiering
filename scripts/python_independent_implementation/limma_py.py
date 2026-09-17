"""Minimal numpy implementations of edgeR TMM normalisation, limma voom, lmFit (weighted) and eBayes.
Follows Robinson & Oshlack 2010 (TMM), Law et al. 2014 (voom) and Smyth 2004 (eBayes)."""
import numpy as np
from scipy import special, stats, optimize


def tmm_factors(counts, logratio_trim=0.3, sum_trim=0.05, a_cutoff=-1e10):
    """counts: genes x samples. Returns normalisation factors (product = 1)."""
    lib = counts.sum(0)
    uq = np.quantile(counts / lib, 0.75, axis=0)
    ref = np.argmin(np.abs(uq - uq.mean()))
    f = np.ones(counts.shape[1])
    for j in range(counts.shape[1]):
        obs, rf = counts[:, j], counts[:, ref]
        nO, nR = lib[j], lib[ref]
        with np.errstate(divide='ignore', invalid='ignore'):
            logR = np.log2((obs / nO) / (rf / nR))
            absE = (np.log2(obs / nO) + np.log2(rf / nR)) / 2
            v = (nO - obs) / nO / obs + (nR - rf) / nR / rf
        ok = np.isfinite(logR) & np.isfinite(absE) & (absE > a_cutoff)
        logR, absE, v = logR[ok], absE[ok], v[ok]
        n = len(logR)
        loL, hiL = np.floor(n * logratio_trim) + 1, n + 1 - np.floor(n * logratio_trim)
        loS, hiS = np.floor(n * sum_trim) + 1, n + 1 - np.floor(n * sum_trim)
        rl = stats.rankdata(logR); rs = stats.rankdata(absE)
        keep = (rl >= loL) & (rl <= hiL) & (rs >= loS) & (rs <= hiS)
        f[j] = 2 ** (np.sum(logR[keep] / v[keep]) / np.sum(1 / v[keep]))
    return f / np.exp(np.mean(np.log(f)))


def _lowess(x, y, frac=0.5, iters=3):
    """Simple LOWESS (tricube, local linear) evaluated at x; returns interpolating function."""
    order = np.argsort(x); xs, ys = x[order], y[order]
    n = len(xs); r = int(np.ceil(frac * n))
    grid = np.linspace(xs.min(), xs.max(), 200)
    rob = np.ones(n)
    for it in range(iters):
        fit = np.empty(len(grid))
        for k, g in enumerate(grid):
            d = np.abs(xs - g); h = np.sort(d)[r - 1] + 1e-12
            w = np.clip(1 - (d / h) ** 3, 0, None) ** 3 * rob
            W = w.sum(); mx = (w * xs).sum() / W; my = (w * ys).sum() / W
            b = (w * (xs - mx) * (ys - my)).sum() / max((w * (xs - mx) ** 2).sum(), 1e-12)
            fit[k] = my + b * (g - mx)
        res = ys - np.interp(xs, grid, fit)
        s = np.median(np.abs(res)) + 1e-12
        rob = np.clip(1 - (res / (6 * s)) ** 2, 0, None) ** 2
    return lambda z: np.interp(z, grid, fit)


def lm_weighted(Y, X, W=None):
    """Y genes x samples, X samples x p, W genes x samples weights. Returns coef, stdev_unscaled, sigma, df."""
    G, n = Y.shape; p = X.shape[1]
    if W is None: W = np.ones_like(Y)
    coef = np.empty((G, p)); su = np.empty((G, p)); sigma = np.empty(G)
    for g in range(G):
        w = W[g]; Xw = X * w[:, None]
        XtWX = X.T @ Xw
        inv = np.linalg.pinv(XtWX)
        b = inv @ (Xw.T @ Y[g])
        r = Y[g] - X @ b
        coef[g] = b; su[g] = np.sqrt(np.diag(inv))
        sigma[g] = np.sqrt((w * r ** 2).sum() / (n - p))
    return coef, su, sigma, n - p


def voom(counts, X, norm_factors=None):
    lib = counts.sum(0) * (norm_factors if norm_factors is not None else 1)
    E = np.log2((counts + 0.5) / (lib + 1) * 1e6)
    coef, su, sigma, df = lm_weighted(E, X)
    Amean = E.mean(1)
    sx = Amean + np.mean(np.log2(lib + 1)) - np.log2(1e6)
    sy = np.sqrt(sigma)
    f = _lowess(sx, sy)
    fitted = coef @ X.T
    fitted_count = fitted + np.log2(lib + 1)[None, :] - np.log2(1e6)
    W = 1 / np.clip(f(fitted_count), 1e-6, None) ** 4
    return E, W


def _trigamma_inverse(x):
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.all(np.abs(dif / y) < 1e-8): break
    return y


def fit_fdist(s2, df):
    ok = np.isfinite(s2) & (s2 > 0)
    z = np.log(s2[ok]); e = z - special.digamma(df / 2) + np.log(df / 2)
    emean = e.mean(); evar = np.mean((e - emean) ** 2) * len(e) / (len(e) - 1) - special.polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s0 = np.exp(emean + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf; s0 = np.exp(emean)
    return d0, s0


def ebayes_test(coef, su, sigma, df, j):
    s2 = sigma ** 2
    d0, s0 = fit_fdist(s2, df)
    if np.isinf(d0): s2p = np.full_like(s2, s0); dft = np.inf
    else: s2p = (d0 * s0 + df * s2) / (d0 + df); dft = d0 + df
    t = coef[:, j] / (su[:, j] * np.sqrt(s2p))
    p = 2 * stats.t.sf(np.abs(t), dft if np.isfinite(dft) else 1e6)
    return t, p, d0, s0


def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1); return out


def filter_by_expr(counts, group_sizes_min, min_count=10, min_total=15):
    lib = counts.sum(0); med = np.median(lib)
    cpm_cut = min_count / med * 1e6
    cpm = counts / lib * 1e6
    return ((cpm >= cpm_cut).sum(1) >= group_sizes_min) & (counts.sum(1) >= min_total)
