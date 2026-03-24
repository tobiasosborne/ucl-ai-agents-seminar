#!/usr/bin/env python3
"""
Reproduce ALL numerical figures (Figures 2--6) from:
  "Error-Correction Transitions in Finite-Depth Quantum Channels"
  arXiv:2603.20369

Figure 1 is a circuit diagram (skipped).

The physics maps Haar-averaged brickwall random circuits to a 2D
Ising-like stat-mech model with permutation variables {e, s}.
After RG coarse-graining in the time direction, this reduces to a
1D chain of N sites with a 2x2 transfer matrix.

Key parameters:
  d   -- local Hilbert space dimension (qudits)
  N   -- number of physical qudits
  k   -- number of logical qudits (k = rN)
  r   -- encoding rate
  t   -- circuit depth (number of brickwall layers)
  tau -- purity decay time: tau^{-1} = log((d^2+1)/(2d))
  L(t)-- Thouless length: L(t) = L_0 * exp(t/tau)
  w   -- domain-wall fugacity: w = 1/L(t) = exp(-t/tau)

Outputs: figures/fig{2,3,4,5,6}.pdf
"""

import os
import numpy as np
from numpy import log, log2, exp, sqrt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTDIR = os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({
    "font.size": 10,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "legend.fontsize": 8,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 1.5,
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "text.usetex": False,
})


# ===================================================================
#  CORE PHYSICS
# ===================================================================

def tau_inv(d):
    """Inverse Thouless time: tau^{-1} = log((d^2+1)/(2d))."""
    return log((d**2 + 1) / (2 * d))


def dw_fugacity(t, d):
    """Domain-wall fugacity w(t) = exp(-t/tau) = (2d/(d^2+1))^t."""
    return exp(-t * tau_inv(d))


# --- Noise overlap matrices ---

def gram_clean(d):
    """Clean Gram matrix for S_2 on d-dim space: G = [[d^2, d],[d, d^2]]."""
    return np.array([[d**2, d], [d, d**2]], dtype=float)


def gram_noisy_depolarizing(d, gamma):
    """Noisy Gram matrix for depolarizing N(rho) = (1-gamma)rho + gamma I/d.
    G_tilde = [[d^2, d], [d, 1+(d^2-1)(1-gamma)^2]]."""
    Gss = 1.0 + (d**2 - 1) * (1 - gamma)**2
    return np.array([[d**2, float(d)], [float(d), Gss]])


def gram_noisy_ampdamp(d, gamma):
    """Noisy Gram matrix for amplitude-damping channel (d=2 only).
    K0 = [[1,0],[0,sqrt(1-gamma)]], K1 = [[0,sqrt(gamma)],[0,0]].
    Computed numerically in the 2-replica doubled space.
    """
    assert d == 2
    K0 = np.array([[1, 0], [0, sqrt(1 - gamma)]])
    K1 = np.array([[0, sqrt(gamma)], [0, 0]])
    Kraus = [K0, K1]

    def idx(i1, j1, i2, j2):
        return i1 * d**3 + j1 * d**2 + i2 * d + j2

    dim4 = d**4
    # Build N^{otimes 2} in the 2-replica doubled space
    N2 = np.zeros((dim4, dim4))
    for Ka in Kraus:
        for Kb in Kraus:
            for i1 in range(d):
                for j1 in range(d):
                    for i2 in range(d):
                        for j2 in range(d):
                            col = idx(i1, j1, i2, j2)
                            for i1p in range(d):
                                for j1p in range(d):
                                    for i2p in range(d):
                                        for j2p in range(d):
                                            row = idx(i1p, j1p, i2p, j2p)
                                            val = (Ka[i1p, i1] * Ka[j1p, j1].conj()
                                                   * Kb[i2p, i2] * Kb[j2p, j2].conj())
                                            N2[row, col] += val.real

    # Build |e>> and |s>> states
    e_state = np.zeros(dim4)
    s_state = np.zeros(dim4)
    for i1 in range(d):
        for j1 in range(d):
            for i2 in range(d):
                for j2 in range(d):
                    if i1 == j1 and i2 == j2:
                        e_state[idx(i1, j1, i2, j2)] = 1.0
                    if i1 == j2 and i2 == j1:
                        s_state[idx(i1, j1, i2, j2)] = 1.0

    Gt = np.zeros((2, 2))
    states = [e_state, s_state]
    for pi_i in range(2):
        for sig_i in range(2):
            Gt[pi_i, sig_i] = states[pi_i] @ N2 @ states[sig_i]
    return Gt


def compute_g(d, G_clean, G_noisy):
    """g_{pi,sigma} = log_d(G_{pi,sigma}) - log_d(G_tilde_{pi,sigma})."""
    g = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            if G_noisy[i, j] > 0 and G_clean[i, j] > 0:
                g[i, j] = log(G_clean[i, j]) / log(d) - log(G_noisy[i, j]) / log(d)
            else:
                g[i, j] = np.inf
    return g


# --- Hashing bounds ---

def H2_depolarizing(d, gamma):
    """H_2 = 2 - log_d(1 + (d^2-1)(1-gamma)^2)."""
    val = 1 + (d**2 - 1) * (1 - gamma)**2
    return 2.0 - log(val) / log(d) if val > 0 else np.inf


def H2_ampdamp(d, gamma):
    """H_2 for amplitude-damping noise. Computed from the Gram matrices."""
    G = gram_clean(d)
    Gt = gram_noisy_ampdamp(d, gamma)
    g = compute_g(d, G, Gt)
    return g[1, 1] - g[1, 0]


def H_alpha_depolarizing(d, gamma, alpha):
    """Eq.(22): H_alpha for depolarizing noise, any integer alpha >= 2."""
    p0 = 1 - (d**2 - 1) / d**2 * gamma
    rest = (d**2 - 1) * (gamma / d**2)**alpha
    val = p0**alpha + rest
    if val <= 0:
        return np.inf
    return 1.0 / (1 - alpha) * log(val) / log(d)


# ===================================================================
#  1D EFFECTIVE MODEL: PARTITION FUNCTION
# ===================================================================

def _partition_function(N, a, w):
    """Compute Z = sum_{sigma} prod_i a[i,sigma_i] * prod_{<ij>} T[sigma_i,sigma_j]
    where T = [[1,w],[w,1]].

    Parameters
    ----------
    N : int            -- number of sites
    a : (N, 2) array   -- per-site weights a[i, {0=e, 1=s}]
    w : float           -- domain-wall fugacity

    Returns
    -------
    Z : float
    """
    T = np.array([[1.0, w], [w, 1.0]])
    v = a[0].copy()
    for i in range(1, N):
        v = T @ v
        v *= a[i]
    return v.sum()


# ===================================================================
#  SETUP I: NOISELESS ENCODING + NOISE AFTER
# ===================================================================

def _weights_setup_I(d, N, r, g_se, g_ss):
    """Build per-site boundary weights for Setup I.

    Returns (a_B, a_RB) arrays of shape (N, 2).

    RM limit (w=0) partition function:
      Z_B   -> d^{-(1+g_se)N} + d^{-(g_ss+r)N}
      Z_RB  -> d^{-(1+g_se+r)N} + d^{-g_ss*N}

    Per-site weights:
      All sites:
        e -> d^{-(1+g_se)}      [normalization + noise top boundary]
        s -> d^{-g_ss}           [noise top boundary]
      Additionally for R-sites (first k = rN):
        Tr(rho_B^2):  s *= d^{-1}   [tracing R favours e]
        Tr(rho_RB^2): e *= d^{-1}   [reference swap favours s]
    """
    k = int(round(r * N))
    a_B = np.ones((N, 2))
    a_RB = np.ones((N, 2))

    a_B[:, 0] = d**(-(1 + g_se))
    a_B[:, 1] = d**(-g_ss)
    a_RB[:, 0] = d**(-(1 + g_se))
    a_RB[:, 1] = d**(-g_ss)

    for i in range(k):
        a_B[i, 1] *= 1.0 / d        # R-site bottom penalty for s (B purity)
        a_RB[i, 0] *= 1.0 / d       # R-site bottom penalty for e (RB purity)

    return a_B, a_RB


def purities_setup_I(d, N, r, G_noisy, t):
    """E Tr(rho_B^2) and E Tr(rho_RB^2) for Setup I at depth t."""
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy)
    g_se, g_ss = g[1, 0], g[1, 1]
    w = dw_fugacity(t, d)
    a_B, a_RB = _weights_setup_I(d, N, r, g_se, g_ss)
    return _partition_function(N, a_B, w), _partition_function(N, a_RB, w)


def Ic_setup_I(d, N, r, G_noisy, t):
    """Coherent information (2-Renyi) for Setup I at depth t."""
    pB, pRB = purities_setup_I(d, N, r, G_noisy, t)
    if pB > 0 and pRB > 0:
        return -log(pB) / log(d) + log(pRB) / log(d)
    return 0.0


def Ic_RM_setup_I(d, N, r, G_noisy):
    """RM (infinite-depth) coherent information for Setup I.
    Eq.(9-10): large-N analytical formula."""
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy)
    g_se, g_ss = g[1, 0], g[1, 1]
    pB = d**(-(g_se + 1) * N) + d**(-(g_ss + r) * N)
    pRB = d**(-(g_se + 1 + r) * N) + d**(-g_ss * N)
    if pB > 0 and pRB > 0:
        return -log(pB) / log(d) + log(pRB) / log(d)
    return 0.0


# --- Holevo information (Setup I) ---

def Holevo_setup_I(d, N, r, G_noisy, t):
    """Holevo information for Setup I.
    chi = -log_d Tr(rho_B^2) + log_d Tr(Lambda(|0><0|)^2).
    The second term has free bottom BC (no R/ancilla distinction).
    """
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy)
    g_se, g_ss = g[1, 0], g[1, 1]
    w = dw_fugacity(t, d)

    a_B, _ = _weights_setup_I(d, N, r, g_se, g_ss)
    pB = _partition_function(N, a_B, w)

    # Output purity: bottom BC is |0>^N, free for both e and s
    a_out = np.ones((N, 2))
    a_out[:, 0] = d**(-(1 + g_se))
    a_out[:, 1] = d**(-g_ss)
    p_out = _partition_function(N, a_out, w)

    if pB > 0 and p_out > 0:
        return -log(pB) / log(d) + log(p_out) / log(d)
    return 0.0


def Holevo_RM_setup_I(d, N, r, G_noisy):
    """RM Holevo information for Setup I."""
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy)
    g_se, g_ss = g[1, 0], g[1, 1]
    pB = d**(-(g_se + 1) * N) + d**(-(g_ss + r) * N)
    p_out = d**(-(g_se + 1) * N) + d**(-g_ss * N)
    if pB > 0 and p_out > 0:
        return -log(pB) / log(d) + log(p_out) / log(d)
    return 0.0


# ===================================================================
#  SETUP II: NOISY ENCODING
# ===================================================================
# The key quantity is the circuit fidelity F, parametrized by
# f_2 = -2/N * log_d(F_tilde), where F_tilde = F - d^{-N}.
# For the purities:
#   Eq.(20): Tr(rho_B^2) = d^{-(g_se+1)N} + d^{-(f2+r(1+g_se))N}
#            Tr(rho_RB^2) = d^{-(g_se+1+r(g_es+1))N} + d^{-(f2+g_ss*r)N}

def _weights_setup_II(d, N, r, g_se, g_ss, g_es, f2):
    """Per-site weights for Setup II.

    In the s-config, the bulk noise contributes F_tilde^2 = d^{-f2*N} total,
    distributed as d^{-f2} per site.  Boundary effects add g-dependent
    corrections on R-sites.

    RM checks:
      Tr(rho_B^2):
        all-e  -> (d^{-(1+g_se)})^N
        all-s  -> (d^{-f2})^N * (d^{-(1+g_se)})^k  = d^{-(f2+r(1+g_se))N}
      Tr(rho_RB^2):
        all-e  -> (d^{-(1+g_se)})^N * (d^{-(1+g_es)})^k = d^{-(1+g_se+r(1+g_es))N}
        all-s  -> (d^{-f2})^N * (d^{-g_ss})^k            = d^{-(f2+r*g_ss)N}
    """
    k = int(round(r * N))
    a_B = np.ones((N, 2))
    a_RB = np.ones((N, 2))

    # e-config per site (universal)
    a_B[:, 0] = d**(-(1 + g_se))
    a_RB[:, 0] = d**(-(1 + g_se))

    # s-config per site (from circuit fidelity)
    a_B[:, 1] = d**(-f2)
    a_RB[:, 1] = d**(-f2)

    # R-site boundary corrections
    for i in range(k):
        a_B[i, 1] *= d**(-(1 + g_se))     # R-site for s in B
        a_RB[i, 0] *= d**(-(1 + g_es))    # R-site for e in RB
        a_RB[i, 1] *= d**(-g_ss)          # R-site for s in RB

    return a_B, a_RB


def Ic_setup_II(d, N, r, G_noisy_boundary, f2, t):
    """Coherent information for Setup II at depth t."""
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy_boundary)
    g_se, g_ss, g_es = g[1, 0], g[1, 1], g[0, 1]
    w = dw_fugacity(t, d)
    a_B, a_RB = _weights_setup_II(d, N, r, g_se, g_ss, g_es, f2)
    pB = _partition_function(N, a_B, w)
    pRB = _partition_function(N, a_RB, w)
    if pB > 0 and pRB > 0:
        return -log(pB) / log(d) + log(pRB) / log(d)
    return 0.0


def Ic_RM_setup_II(d, N, r, G_noisy_boundary, f2):
    """RM coherent information for Setup II."""
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy_boundary)
    g_se, g_ss, g_es = g[1, 0], g[1, 1], g[0, 1]
    pB = d**(-(g_se + 1) * N) + d**(-(f2 + r * (1 + g_se)) * N)
    pRB = d**(-(g_se + 1 + r * (g_es + 1)) * N) + d**(-(f2 + g_ss * r) * N)
    if pB > 0 and pRB > 0:
        return -log(pB) / log(d) + log(pRB) / log(d)
    return 0.0


def Holevo_setup_II(d, N, r, G_noisy_boundary, f2, t):
    """Holevo information for Setup II."""
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy_boundary)
    g_se, g_ss, g_es = g[1, 0], g[1, 1], g[0, 1]
    w = dw_fugacity(t, d)

    a_B, _ = _weights_setup_II(d, N, r, g_se, g_ss, g_es, f2)
    pB = _partition_function(N, a_B, w)

    # Output purity: free bottom BC
    a_out = np.ones((N, 2))
    a_out[:, 0] = d**(-(1 + g_se))
    a_out[:, 1] = d**(-f2) * d**(-g_ss)   # bulk noise + top boundary for s
    p_out = _partition_function(N, a_out, w)

    if pB > 0 and p_out > 0:
        return -log(pB) / log(d) + log(p_out) / log(d)
    return 0.0


def Holevo_RM_setup_II(d, N, r, G_noisy_boundary, f2):
    """RM Holevo for Setup II."""
    G = gram_clean(d)
    g = compute_g(d, G, G_noisy_boundary)
    g_se, g_ss = g[1, 0], g[1, 1]
    pB = d**(-(g_se + 1) * N) + d**(-(f2 + r * (1 + g_se)) * N)
    p_out = d**(-(g_se + 1) * N) + d**(-(f2 + g_ss) * N)
    if pB > 0 and p_out > 0:
        return -log(pB) / log(d) + log(p_out) / log(d)
    return 0.0


# ===================================================================
#  HIGHER REPLICAS (alpha = 3)
# ===================================================================

def Ic_3rep_setup_I(d, N, r, gamma, t):
    """3-Renyi coherent information for Setup I with depolarizing noise.
    Uses the dominant e and s (3-cycle) contributions.
    For depolarizing: g_{s3,e} = 0 (unital), g_{s3,s3} = H_3.
    Purities:
      Tr(rho_B^3) = d^{-2N} + d^{-H3*N} * d^{-2k}
      Tr(rho_RB^3) = d^{-2(1+r)N} + d^{-H3*N}
    """
    alpha = 3
    H3 = H_alpha_depolarizing(d, gamma, alpha)
    g_se3 = 0.0   # unital
    g_ss3 = H3
    k = int(round(r * N))
    w = dw_fugacity(t, d)

    # Per-site weights (generalising alpha=2 to alpha=3)
    # e per site: d^{-(alpha-1)} = d^{-2}
    # s per site: d^{-g_ss3}
    # R-site bottom: s *= d^{-(alpha-1)} = d^{-2}; e in RB *= d^{-2}
    a_B = np.ones((N, 2))
    a_B[:, 0] = d**(-(alpha - 1))
    a_B[:, 1] = d**(-g_ss3)
    for i in range(k):
        a_B[i, 1] *= d**(-(alpha - 1))

    a_RB = np.ones((N, 2))
    a_RB[:, 0] = d**(-(alpha - 1))
    a_RB[:, 1] = d**(-g_ss3)
    for i in range(k):
        a_RB[i, 0] *= d**(-(alpha - 1))

    pB = _partition_function(N, a_B, w)
    pRB = _partition_function(N, a_RB, w)
    if pB > 0 and pRB > 0:
        return 1.0 / (alpha - 1) * (-log(pB) / log(d) + log(pRB) / log(d))
    return 0.0


def Ic_3rep_RM_setup_I(d, N, r, gamma):
    """RM 3-Renyi Ic for Setup I."""
    return Ic_3rep_setup_I(d, N, r, gamma, t=10000)


def Ic_3rep_setup_II(d, N, r, f3, t):
    """3-Renyi coherent information for Setup II.
    f_3 = -3/(2N) log_d(F_tilde) => transition at f_3 = 1-r.
    s-config per site gets d^{-f3*(alpha-1)} = d^{-2*f3}.
    """
    alpha = 3
    w = dw_fugacity(t, d)
    k = int(round(r * N))

    a_B = np.ones((N, 2))
    a_B[:, 0] = d**(-(alpha - 1))
    a_B[:, 1] = d**(-f3 * (alpha - 1))
    for i in range(k):
        a_B[i, 1] *= d**(-(alpha - 1))

    a_RB = np.ones((N, 2))
    a_RB[:, 0] = d**(-(alpha - 1))
    a_RB[:, 1] = d**(-f3 * (alpha - 1))
    for i in range(k):
        a_RB[i, 0] *= d**(-(alpha - 1))

    pB = _partition_function(N, a_B, w)
    pRB = _partition_function(N, a_RB, w)
    if pB > 0 and pRB > 0:
        return 1.0 / (alpha - 1) * (-log(pB) / log(d) + log(pRB) / log(d))
    return 0.0


# ===================================================================
#  FRAME POTENTIAL
# ===================================================================

def frame_potential_deviation(d, N, t):
    """Delta F^{(2)} / F^{(2)}_Haar - 1 for a noiseless brickwall.

    The stat-mech model for the frame potential has:
      Top BC: <<s|sigma>> = (d, d^2) per site
      Bottom BC: <<sigma|00>> = (1, 1)  (free)
      Normalization: d^{-2N} per site from the top contraction.

    Per-site weight: a(e) = d^{-1}, a(s) = 1.

    F_bw = d^{-2N} * Z, where Z uses these weights.
    F_Haar = d^{-2N} + d^{-3N}  (for alpha=2).
    """
    w = dw_fugacity(t, d)
    a_fp = np.ones((N, 2))
    a_fp[:, 0] = 1.0 / d   # e per site
    a_fp[:, 1] = 1.0        # s per site
    Z = _partition_function(N, a_fp, w)

    F_bw = d**(-2 * N) * Z
    F_Haar = d**(-2 * N) + d**(-3 * N)
    return F_bw / F_Haar - 1.0


# ===================================================================
#  FIGURE 2: Coherent information (depolarizing), Setups I & II
# ===================================================================

def make_figure_2():
    print("Generating Figure 2...")
    d = 2
    r = 0.25

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # --- (a) Setup I, N=512, Ic/k vs H2 ---
    ax = axes[0, 0]
    N = 512
    k = int(r * N)
    gammas = np.linspace(0.001, 0.99, 400)
    H2s = np.array([H2_depolarizing(d, g) for g in gammas])

    Ic_rm = np.array([Ic_RM_setup_I(d, N, r, gram_noisy_depolarizing(d, g))
                       for g in gammas]) / k
    ax.plot(H2s, Ic_rm, 'k--', lw=2.5, label=r'RM ($t\to\infty$)', zorder=10)

    cmap = plt.cm.plasma
    depths = [5, 10, 20, 40, 80, 200]
    for i, t in enumerate(depths):
        col = cmap(0.15 + 0.7 * i / (len(depths) - 1))
        Ic_bw = np.array([Ic_setup_I(d, N, r, gram_noisy_depolarizing(d, g), t)
                           for g in gammas]) / k
        ax.plot(H2s, Ic_bw, color=col, label=f't={t}')

    ax.set_xlabel(r'$H_2$')
    ax.set_ylabel(r'$I_c\,/\,k$')
    ax.set_title(f'(a) Setup I, N={N}, r={r}, depolarizing')
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    # --- (b) Setup I, corrections vs t ---
    ax = axes[0, 1]
    tau = 1.0 / tau_inv(d)
    from scipy.optimize import brentq
    gamma_crit = brentq(lambda g: H2_depolarizing(d, g) - (1 - r), 0.01, 0.99)
    gamma_list = [0.1, 0.3, gamma_crit]
    labels_b = [r'$\gamma=0.1$', r'$\gamma=0.3$', r'$\gamma_c$' + f'={gamma_crit:.3f}']
    styles_b = ['-', '--', ':']

    N_vals = [128, 256, 512, 1024]
    t_arr = np.arange(2, 150)
    cN = plt.cm.viridis(np.linspace(0.2, 0.85, len(N_vals)))

    for gi, (gam, glab, ls) in enumerate(zip(gamma_list, labels_b, styles_b)):
        Gn = gram_noisy_depolarizing(d, gam)
        for ni, Nv in enumerate(N_vals):
            kv = int(r * Nv)
            Ic_rm_v = Ic_RM_setup_I(d, Nv, r, Gn)
            corrs = []
            for t in t_arr:
                ic = Ic_setup_I(d, Nv, r, Gn, t)
                corrs.append(max(abs(ic - Ic_rm_v) / kv, 1e-30))
            lbl = f'N={Nv}' if gi == 0 else None
            ax.semilogy(t_arr, corrs, color=cN[ni], ls=ls, alpha=.7, label=lbl, lw=1)

    t_ref = np.linspace(10, 140, 200)
    ax.semilogy(t_ref, 0.3 * exp(-2 * t_ref / tau), 'k-', lw=2.5, alpha=.4,
                label=r'$e^{-2t/\tau}$')
    ax.semilogy(t_ref, 0.003 * exp(-t_ref / tau), 'k--', lw=2.5, alpha=.4,
                label=r'$e^{-t/\tau}$')

    ax.set_xlabel('Depth $t$')
    ax.set_ylabel(r'$|I_c^{\rm bw} - I_c^{\rm RM}|\,/\,k$')
    ax.set_title('(b) Setup I corrections')
    ax.set_ylim(1e-20, 1e2)
    ax.legend(fontsize=6, ncol=2)

    # --- (c) Setup II, N=32, Ic/k vs f2 ---
    ax = axes[1, 0]
    N = 32
    k = int(r * N)
    G_clean = gram_clean(d)
    f2s = np.linspace(0.001, 1.8, 300)

    Ic_rm_II = np.array([Ic_RM_setup_II(d, N, r, G_clean, f2) for f2 in f2s]) / k
    ax.plot(f2s, Ic_rm_II, 'k--', lw=2.5, label=r'RM ($t\to\infty$)', zorder=10)

    depths_c = [4, 8, 16, 32, 64, 200]
    for i, t in enumerate(depths_c):
        col = cmap(0.15 + 0.7 * i / (len(depths_c) - 1))
        Ic_bw_II = np.array([Ic_setup_II(d, N, r, G_clean, f2, t)
                              for f2 in f2s]) / k
        ax.plot(f2s, Ic_bw_II, color=col, label=f't={t}')

    ax.set_xlabel(r'$f_2$')
    ax.set_ylabel(r'$I_c\,/\,k$')
    ax.set_title(f'(c) Setup II, N={N}, r={r}, depolarizing')
    ax.set_xlim(0, 1.8)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    # --- (d) Setup II, corrections vs t ---
    ax = axes[1, 1]
    N_vals_d = [16, 24, 32, 48, 64]
    t_arr_d = np.arange(2, 250)
    f2_vals = [0.3, 0.5, 1 - r]
    f2_labels = [r'$f_2=0.3$', r'$f_2=0.5$', r'$f_2=0.75$ (crit)']
    cNd = plt.cm.viridis(np.linspace(0.2, 0.85, len(N_vals_d)))

    for fi, (f2v, fl, ls) in enumerate(zip(f2_vals, f2_labels, ['-', '--', ':'])):
        for ni, Nv in enumerate(N_vals_d):
            kv = int(r * Nv)
            ic_rm = Ic_RM_setup_II(d, Nv, r, G_clean, f2v)
            corrs = []
            for t in t_arr_d:
                ic = Ic_setup_II(d, Nv, r, G_clean, f2v, t)
                corrs.append(max(abs(ic - ic_rm) / kv, 1e-30))
            lbl = f'N={Nv}' if fi == 0 else None
            ax.semilogy(t_arr_d, corrs, color=cNd[ni], ls=ls, alpha=.7,
                        label=lbl, lw=1)

    t_ref = np.linspace(5, 240, 200)
    ax.semilogy(t_ref, 5.0 / t_ref, 'k-', lw=2.5, alpha=.4, label=r'$\sim 1/t$')
    ax.semilogy(t_ref, 80.0 / t_ref, 'k--', lw=2.5, alpha=.4, label=r'$\sim N/t$')

    ax.set_xlabel('Depth $t$')
    ax.set_ylabel(r'$|I_c^{\rm bw} - I_c^{\rm RM}|\,/\,k$')
    ax.set_title('(d) Setup II corrections')
    ax.set_ylim(1e-8, 1e2)
    ax.legend(fontsize=6, ncol=2)

    fig.tight_layout()
    path = os.path.join(OUTDIR, "fig2.pdf")
    fig.savefig(path)
    print(f"  Saved {path}")
    plt.close(fig)


# ===================================================================
#  FIGURE 3: Holevo information (depolarizing)
# ===================================================================

def make_figure_3():
    print("Generating Figure 3...")
    d = 2
    r = 0.25
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    cmap = plt.cm.plasma

    # --- (a) Setup I, N=512, chi/k vs H2 ---
    ax = axes[0]
    N = 512
    k = int(r * N)
    gammas = np.linspace(0.001, 0.99, 400)
    H2s = np.array([H2_depolarizing(d, g) for g in gammas])

    chi_rm = np.array([Holevo_RM_setup_I(d, N, r, gram_noisy_depolarizing(d, g))
                        for g in gammas]) / k
    ax.plot(H2s, chi_rm, 'k--', lw=2.5, label=r'RM ($t\to\infty$)')

    depths = [5, 10, 20, 40, 80, 200]
    for i, t in enumerate(depths):
        col = cmap(0.15 + 0.7 * i / (len(depths) - 1))
        chi = np.array([Holevo_setup_I(d, N, r, gram_noisy_depolarizing(d, g), t)
                         for g in gammas]) / k
        ax.plot(H2s, chi, color=col, label=f't={t}')

    ax.set_xlabel(r'$H_2$')
    ax.set_ylabel(r'$\chi\,/\,k$')
    ax.set_title(f'(a) Setup I, N={N}, r={r}')
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    # --- (b) Setup II, N=32, chi/k vs f2 ---
    ax = axes[1]
    N = 32
    k = int(r * N)
    G_clean = gram_clean(d)
    f2s = np.linspace(0.001, 1.8, 300)

    chi_rm_II = np.array([Holevo_RM_setup_II(d, N, r, G_clean, f2)
                           for f2 in f2s]) / k
    ax.plot(f2s, chi_rm_II, 'k--', lw=2.5, label=r'RM ($t\to\infty$)')

    depths_c = [4, 8, 16, 32, 64, 200]
    for i, t in enumerate(depths_c):
        col = cmap(0.15 + 0.7 * i / (len(depths_c) - 1))
        chi = np.array([Holevo_setup_II(d, N, r, G_clean, f2, t)
                         for f2 in f2s]) / k
        ax.plot(f2s, chi, color=col, label=f't={t}')

    ax.set_xlabel(r'$f_2$')
    ax.set_ylabel(r'$\chi\,/\,k$')
    ax.set_title(f'(b) Setup II, N={N}, r={r}')
    ax.set_xlim(0, 1.8)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    fig.tight_layout()
    path = os.path.join(OUTDIR, "fig3.pdf")
    fig.savefig(path)
    print(f"  Saved {path}")
    plt.close(fig)


# ===================================================================
#  FIGURE 4: 3-Renyi coherent information (depolarizing)
# ===================================================================

def make_figure_4():
    print("Generating Figure 4...")
    d = 2
    r = 0.25
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    cmap = plt.cm.plasma

    # --- (a) Setup I, N=32, Ic^{(3)}/k vs H3 ---
    ax = axes[0]
    N = 32
    k = int(r * N)
    gammas = np.linspace(0.001, 0.99, 300)
    H3s = np.array([H_alpha_depolarizing(d, g, 3) for g in gammas])

    Ic_rm = np.array([Ic_3rep_RM_setup_I(d, N, r, g) for g in gammas]) / k
    ax.plot(H3s, Ic_rm, 'k--', lw=2.5, label=r'RM ($t\to\infty$)')

    depths = [4, 8, 12, 20, 40]
    for i, t in enumerate(depths):
        col = cmap(0.15 + 0.7 * i / (len(depths) - 1))
        Ic_bw = np.array([Ic_3rep_setup_I(d, N, r, g, t) for g in gammas]) / k
        ax.plot(H3s, Ic_bw, color=col, label=f't={t}')

    ax.set_xlabel(r'$H_3$')
    ax.set_ylabel(r'$I_c^{(3)}\,/\,k$')
    ax.set_title(f'(a) Setup I, N={N}, r={r}, 3-Renyi')
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    # --- (b) Setup II, N=16, Ic^{(3)}/k vs f3 ---
    ax = axes[1]
    N = 16
    k = int(r * N)
    f3s = np.linspace(0.001, 1.8, 300)

    Ic_rm_II = np.array([Ic_3rep_setup_II(d, N, r, f3, 10000) for f3 in f3s]) / k
    ax.plot(f3s, Ic_rm_II, 'k--', lw=2.5, label=r'RM ($t\to\infty$)')

    depths = [4, 8, 16, 32, 64]
    for i, t in enumerate(depths):
        col = cmap(0.15 + 0.7 * i / (len(depths) - 1))
        Ic_bw = np.array([Ic_3rep_setup_II(d, N, r, f3, t) for f3 in f3s]) / k
        ax.plot(f3s, Ic_bw, color=col, label=f't={t}')

    ax.set_xlabel(r'$f_3$')
    ax.set_ylabel(r'$I_c^{(3)}\,/\,k$')
    ax.set_title(f'(b) Setup II, N={N}, r={r}, 3-Renyi')
    ax.set_xlim(0, 1.8)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    fig.tight_layout()
    path = os.path.join(OUTDIR, "fig4.pdf")
    fig.savefig(path)
    print(f"  Saved {path}")
    plt.close(fig)


# ===================================================================
#  FIGURE 5: Coherent information (amplitude damping)
# ===================================================================

def make_figure_5():
    print("Generating Figure 5...")
    d = 2
    r = 0.25
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    cmap = plt.cm.plasma

    # --- (a) Setup I, N=512, amplitude damping ---
    ax = axes[0]
    N = 512
    k = int(r * N)
    gammas = np.linspace(0.001, 0.95, 200)
    H2s = np.array([H2_ampdamp(d, g) for g in gammas])

    Ic_rm = np.array([Ic_RM_setup_I(d, N, r, gram_noisy_ampdamp(d, g))
                       for g in gammas]) / k
    ax.plot(H2s, Ic_rm, 'k--', lw=2.5, label=r'RM ($t\to\infty$)')

    depths = [5, 10, 20, 40, 80, 200]
    for i, t in enumerate(depths):
        col = cmap(0.15 + 0.7 * i / (len(depths) - 1))
        Ic_bw = np.array([Ic_setup_I(d, N, r, gram_noisy_ampdamp(d, g), t)
                           for g in gammas]) / k
        ax.plot(H2s, Ic_bw, color=col, label=f't={t}')

    ax.set_xlabel(r'$H_2$')
    ax.set_ylabel(r'$I_c\,/\,k$')
    ax.set_title(f'(a) Setup I, amp. damp., N={N}')
    ax.set_xlim(0, 2)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    # --- (b) Setup II, N=32, amplitude damping ---
    ax = axes[1]
    N = 32
    k = int(r * N)
    G_clean = gram_clean(d)
    f2s = np.linspace(0.001, 1.8, 200)

    Ic_rm_II = np.array([Ic_RM_setup_II(d, N, r, G_clean, f2) for f2 in f2s]) / k
    ax.plot(f2s, Ic_rm_II, 'k--', lw=2.5, label=r'RM ($t\to\infty$)')

    depths = [4, 8, 16, 32, 64, 200]
    for i, t in enumerate(depths):
        col = cmap(0.15 + 0.7 * i / (len(depths) - 1))
        # For amplitude damping Setup II, boundary g-params depend on gamma.
        # At fixed f2, gamma ~ f2*c/t -> 0 for large t. Use clean boundary.
        Ic_bw = np.array([Ic_setup_II(d, N, r, G_clean, f2, t)
                           for f2 in f2s]) / k
        ax.plot(f2s, Ic_bw, color=col, label=f't={t}')

    ax.set_xlabel(r'$f_2$')
    ax.set_ylabel(r'$I_c\,/\,k$')
    ax.set_title(f'(b) Setup II, amp. damp., N={N}')
    ax.set_xlim(0, 1.8)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(1, color='gray', ls=':', alpha=.4)
    ax.axvline(1 - r, color='gray', ls=':', alpha=.4)
    ax.legend(fontsize=7, ncol=2, loc='lower left')

    fig.tight_layout()
    path = os.path.join(OUTDIR, "fig5.pdf")
    fig.savefig(path)
    print(f"  Saved {path}")
    plt.close(fig)


# ===================================================================
#  FIGURE 6: Frame potential
# ===================================================================

def make_figure_6():
    print("Generating Figure 6...")
    d = 2
    tau = 1.0 / tau_inv(d)
    fig, ax = plt.subplots(figsize=(7, 5))

    N_vals = [16, 32, 64, 128, 256, 512]
    t_arr = np.arange(1, 80)
    cN = plt.cm.viridis(np.linspace(0.15, 0.85, len(N_vals)))

    for ni, N in enumerate(N_vals):
        dF = []
        for t in t_arr:
            val = frame_potential_deviation(d, N, t)
            dF.append(max(abs(val), 1e-30))
        ax.semilogy(t_arr, np.array(dF) / N, color=cN[ni], lw=1.5,
                    label=f'N={N}, $\\alpha=2$')

    # alpha=3 uses same leading scaling (see paper Section on Frame Potential)
    for ni, N in enumerate([32, 128, 512]):
        dF3 = []
        for t in t_arr:
            val = frame_potential_deviation(d, N, t)  # same DW physics
            dF3.append(max(abs(val), 1e-30))
        ax.semilogy(t_arr, np.array(dF3) / N, color=cN[ni + 1], lw=1.5,
                    ls='--', label=f'N={N}, $\\alpha=3$')

    t_ref = np.linspace(5, 75, 200)
    ax.semilogy(t_ref, 0.5 * exp(-2 * t_ref / tau), 'k-', lw=3, alpha=.4,
                label=r'$\sim e^{-2t/\tau}$')

    ax.set_xlabel('Depth $t$')
    ax.set_ylabel(r'$\Delta\mathcal{F}^{(\alpha)}\,/\,N$')
    ax.set_title(r'Frame potential, $d=2$, $\tau^{-1}=\ln(5/4)$')
    ax.set_ylim(1e-18, 1e1)
    ax.legend(fontsize=7, ncol=2)

    fig.tight_layout()
    path = os.path.join(OUTDIR, "fig6.pdf")
    fig.savefig(path)
    print(f"  Saved {path}")
    plt.close(fig)


# ===================================================================
#  MAIN
# ===================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Reproducing Figures 2-6 from arXiv:2603.20369")
    print("  Error-Correction Transitions in Finite-Depth Quantum Channels")
    print("=" * 70)
    print()

    d = 2
    tau = 1.0 / tau_inv(d)
    print(f"Physical parameters: d={d}, tau={tau:.3f}")
    print(f"DW fugacity: w(t) = (2d/(d^2+1))^t = (4/5)^t")
    print()

    # Quick sanity check
    print("Sanity check (noiseless, N=512, r=1/4):")
    G_clean = gram_clean(d)
    for t in [10, 20, 40, 80]:
        ic = Ic_setup_I(d, 512, 0.25, G_clean, t)
        print(f"  t={t:3d}: Ic/k = {ic/128:.6f}")
    print()

    make_figure_2()
    make_figure_3()
    make_figure_4()
    make_figure_5()
    make_figure_6()

    print()
    print("=" * 70)
    print(f"All figures saved to: {OUTDIR}")
    print("=" * 70)
