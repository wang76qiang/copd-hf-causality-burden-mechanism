#!/usr/bin/env python3
"""MR engine + harmonisation utilities (pure Python, no external MR deps).
Validated against published discovery IVW OR 1.151 (95% CI 1.083-1.223), I2=52.3%.
"""
import math, csv, json
from scipy import stats, optimize
import numpy as np

def comp(a):
    return a.translate(str.maketrans("ACGT", "TGCA"))

def wald_ivw(bx, by, sey):
    """Wald ratios + fixed/random-effect IVW + Cochran Q/I2."""
    wald = [y/x for x, y in zip(bx, by)]
    se_w = [s/abs(x) for s, x in zip(sey, bx)]
    w = np.array([1/s**2 for s in se_w])
    b = np.array(wald)
    est_f = float(np.sum(w*b)/np.sum(w))
    se_f = math.sqrt(1/np.sum(w))
    Q = float(np.sum(w*(b-est_f)**2))
    df = len(b)-1
    Q_p = float(stats.chi2.sf(Q, df))
    I2 = max(0.0, (Q-df)/Q) if Q > 0 else 0.0
    tau2 = max(0.0, (Q-df)/(np.sum(w) - np.sum(w**2)/np.sum(w)))
    w2 = 1/(np.array(se_w)**2 + tau2)
    est_r = float(np.sum(w2*b)/np.sum(w2))
    se_r = math.sqrt(1/np.sum(w2))
    return {
        "n": len(b),
        "ivw_fe": est_f, "ivw_fe_se": se_f, "ivw_fe_p": float(2*stats.norm.sf(abs(est_f/se_f))),
        "ivw_re": est_r, "ivw_re_se": se_r, "ivw_re_p": float(2*stats.norm.sf(abs(est_r/se_r))),
        "Q": Q, "Q_df": df, "Q_p": Q_p, "I2": I2, "tau2": tau2,
        "wald": wald, "wald_se": se_w,
    }

def mr_egger(bx, by, sey):
    """MR-Egger with exposure betas oriented positive."""
    sgn = np.sign(bx)
    x = np.abs(np.array(bx)); y = np.array(by)*sgn
    w = 1/np.array(sey)**2
    X = np.column_stack([np.ones_like(x), x])
    W = np.diag(w)
    XtW = X.T @ W
    beta = np.linalg.solve(XtW @ X, XtW @ y)
    resid = y - X @ beta
    dfree = len(x) - 2
    sigma2 = max(1.0, float((resid * w * resid).sum() / dfree))  # multiplicative overdispersion, min 1
    cov = sigma2 * np.linalg.inv(XtW @ X)
    se_int, se_slo = math.sqrt(cov[0,0]), math.sqrt(cov[1,1])
    return {
        "intercept": float(beta[0]), "intercept_se": se_int,
        "intercept_p": float(2*stats.t.sf(abs(beta[0]/se_int), dfree)),
        "slope": float(beta[1]), "slope_se": se_slo,
        "slope_p": float(2*stats.t.sf(abs(beta[1]/se_slo), dfree)),
    }

def weighted_median(bx, by, sey, nboot=5000, seed=42):
    wald = np.array([y/x for x, y in zip(bx, by)])
    se_w = np.array([s/abs(x) for s, x in zip(sey, bx)])
    w = 1/se_w**2
    order = np.argsort(wald)
    wald_s, w_s = wald[order], w[order]
    cum = np.cumsum(w_s) - 0.5*w_s
    cum /= np.sum(w_s)
    est = float(np.interp(0.5, cum, wald_s))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(nboot):
        bsamp = rng.normal(wald, se_w)
        o = np.argsort(bsamp); bs, ws = bsamp[o], w[o]
        c = np.cumsum(ws) - 0.5*ws; c /= np.sum(ws)
        boots.append(np.interp(0.5, c, bs))
    se = float(np.std(boots))
    return {"wm": est, "wm_se": se, "wm_p": float(2*stats.norm.sf(abs(est/se)))}

def full_mr(bx, by, sey):
    out = wald_ivw(bx, by, sey)
    out.update(mr_egger(bx, by, sey))
    out.update(weighted_median(bx, by, sey))
    out.pop("wald"); out.pop("wald_se")
    return out

def winners_curse_conditional(b, se, z_thresh=5.45131):
    """Zhong & Prentice style conditional-likelihood correction for selection |Z|>z_thresh.
    MLE of true effect mu given observed b ~ N(mu, se^2) conditional on selection."""
    def nll(mu):
        # conditional density: phi((b-mu)/se)/se / P(|N(mu,se)|/se > z)
        ll = stats.norm.logpdf(b, mu, se)
        tail = stats.norm.cdf(-z_thresh*se - mu, 0, se) + stats.norm.sf(z_thresh*se - mu, 0, se)
        tail = max(tail, 1e-300)
        return -(ll - math.log(tail))
    r = optimize.minimize_scalar(nll, bounds=(b - 10*se, b + 10*se), method="bounded")
    return float(r.x)

def harmonise_outcome(ea_x, oa_x, eaf_x, ea_y, oa_y, beta_y, af_y=None):
    """Align outcome beta to exposure effect allele. Returns (beta_aligned, action) or (None, reason)."""
    if ea_y == ea_x and oa_y == oa_x:
        b, act = beta_y, "aligned"
        ea_freq = af_y
    elif ea_y == oa_x and oa_y == ea_x:
        b, act = -beta_y, "flipped"
        ea_freq = (1 - af_y) if af_y is not None else None
    elif comp(ea_y) == ea_x and comp(oa_y) == oa_x:
        b, act = beta_y, "strand_aligned"
        ea_freq = af_y
    elif comp(ea_y) == oa_x and comp(oa_y) == ea_x:
        b, act = -beta_y, "strand_flipped"
        ea_freq = (1 - af_y) if af_y is not None else None
    else:
        return None, "allele_mismatch"
    pal = frozenset([ea_x, oa_x]) in (frozenset("AT"), frozenset("CG"))
    if pal and eaf_x is not None and ea_freq is not None:
        if abs(ea_freq - eaf_x) > 0.15:
            return None, "ambiguous_palindromic"
    return b, act + (";pal" if pal else "")
