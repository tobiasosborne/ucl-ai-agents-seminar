#!/usr/bin/env python3
"""Reproduce Figure 2 using EXACT 2D stat-mech computation."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from exact_statmech import (coherent_info, Ic_RM, H2_depolarizing,
                             tau_brickwall, purity_exact)
import time

d = 2
r = 0.25
tau = tau_brickwall(d)

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# ─────────────────────────────────────────────────────
# Panel (a): Setup I, Ic/k vs H2
# ─────────────────────────────────────────────────────
ax = axes[0, 0]
N = 16  # exact feasible (paper uses 512; qualitative match expected)
k = int(r * N)
gammas = np.linspace(0, 0.45, 60)
H2s = np.array([H2_depolarizing(d, g) for g in gammas])

depths = [2, 4, 6, 8, 10, 14, 18, 22, 26, 30]
cmap_norm = Normalize(vmin=min(depths), vmax=max(depths))
cmap = plt.cm.cool

# RM curve
Ic_rm_vals = np.array([Ic_RM(d, N, r, g) / k for g in gammas])
ax.plot(H2s, Ic_rm_vals, 'k--', lw=2.5, label=r'$I_c^{\rm RM}$', zorder=10)

t0 = time.time()
for t_val in depths:
    Ic_vals = []
    for g in gammas:
        Ic_vals.append(coherent_info(d, N, r, t_val, g, 'depolarizing', 'I') / k)
    ax.plot(H2s, Ic_vals, color=cmap(cmap_norm(t_val)), lw=1.3)
elapsed = time.time() - t0
print(f"Panel (a): {elapsed:.1f}s for N={N}")

sm = ScalarMappable(cmap=cmap, norm=cmap_norm)
sm.set_array([])
cb = fig.colorbar(sm, ax=ax, label='Depth $t$', shrink=0.8)

ax.axhline(0, color='gray', ls=':', lw=0.5)
ax.axvline(1 - r, color='gray', ls=':', lw=0.5, alpha=0.5)
ax.set_xlabel(r'$H_2$')
ax.set_ylabel(r'$I_c^{\rm bw}/k$')
ax.set_xlim(0, 1.6)
ax.set_ylim(-1, 1.05)
ax.set_title(f'(a) Setup I, $N={N}$, $r={r}$')
ax.legend(loc='upper right', fontsize=9)

# ─────────────────────────────────────────────────────
# Panel (b): Setup I, corrections vs t
# ─────────────────────────────────────────────────────
ax = axes[0, 1]
Ns = [8, 12, 16, 20]
colors = ['C0', 'C1', 'C2', 'C3']
gammas_corr = [0.1, 0.25]  # in recovery regime
markers = ['o', 's']

for gi, gamma in enumerate(gammas_corr):
    for ni, N_val in enumerate(Ns):
        k_val = int(r * N_val)
        Ic_rm_val = Ic_RM(d, N_val, r, gamma)
        ts = np.arange(2, 28, 2)
        corrs = []
        for t_val in ts:
            Ic_bw = coherent_info(d, N_val, r, t_val, gamma, 'depolarizing', 'I')
            corrs.append(abs(Ic_bw - Ic_rm_val) / k_val)
        corrs = np.array(corrs)
        mask = corrs > 1e-15
        label = f'N={N_val}' if gi == 0 else ''
        ax.semilogy(ts[mask], corrs[mask], markers[gi] + '-',
                   color=colors[ni], ms=4, lw=1, label=label)

# Reference lines
t_ref = np.linspace(4, 26, 100)
ax.semilogy(t_ref, 0.3 * np.exp(-t_ref / tau), 'k:', lw=1.5, alpha=0.6, label=r'$e^{-t/\tau}$')
ax.semilogy(t_ref, 0.3 * np.exp(-2 * t_ref / tau), 'k-.', lw=1.5, alpha=0.6, label=r'$e^{-2t/\tau}$')

ax.set_xlabel('$t$')
ax.set_ylabel(r'$|I_c^{\rm bw} - I_c^{\rm RM}|/k$')
ax.set_title('(b) Setup I corrections')
ax.legend(fontsize=8, ncol=2)
print(f"Panel (b) done")

# ─────────────────────────────────────────────────────
# Panel (c): Setup II, Ic/k vs f2
# ─────────────────────────────────────────────────────
ax = axes[1, 0]
N = 12  # smaller for Setup II (more layers needed)
k = int(r * N)
depths_II = [2, 4, 6, 8, 12, 16, 20]
gammas_II = np.linspace(0.01, 0.4, 40)
cmap_norm2 = Normalize(vmin=min(depths_II), vmax=max(depths_II))

t0 = time.time()
for t_val in depths_II:
    f2_vals = []
    Ic_vals = []
    for g in gammas_II:
        # f2 from circuit fidelity (approximate: F ~ ratio^{Nt/2})
        ratio_ss = (1 + (d**2 - 1) * (1 - g)**2) / d**2
        log_Ft = (N * t_val / 2) * np.log(max(ratio_ss, 1e-300))
        f2 = -2 / N * log_Ft / np.log(d)
        f2_vals.append(f2)

        Ic_bw = coherent_info(d, N, r, t_val, g, 'depolarizing', 'II')
        Ic_vals.append(Ic_bw / k)

    ax.plot(f2_vals, Ic_vals, color=cmap(cmap_norm2(t_val)), lw=1.3)

elapsed = time.time() - t0
print(f"Panel (c): {elapsed:.1f}s for N={N}")

ax.axhline(0, color='gray', ls=':', lw=0.5)
ax.axvline(1 - r, color='gray', ls=':', lw=0.5, alpha=0.5)
ax.set_xlabel(r'$f_2$')
ax.set_ylabel(r'$I_c^{\rm bw}/k$')
ax.set_xlim(0, 1.5)
ax.set_ylim(-1, 1.05)
ax.set_title(f'(c) Setup II, $N={N}$, $r={r}$')

sm2 = ScalarMappable(cmap=cmap, norm=cmap_norm2)
sm2.set_array([])
fig.colorbar(sm2, ax=ax, label='Depth $t$', shrink=0.8)

# ─────────────────────────────────────────────────────
# Panel (d): Setup II, corrections vs t
# ─────────────────────────────────────────────────────
ax = axes[1, 1]
Ns_II = [8, 12, 16]
colors_II = ['C0', 'C1', 'C2']

for ni, N_val in enumerate(Ns_II):
    k_val = int(r * N_val)
    # Fixed f2 ≈ 0.5 by adjusting gamma: gamma ≈ f2/(t * const)
    ts = np.arange(4, 30, 2)
    corrs = []
    for t_val in ts:
        gamma = 0.5 / t_val  # scale gamma to keep f2 roughly fixed
        Ic_rm_val = Ic_RM(d, N_val, r, gamma)
        Ic_bw = coherent_info(d, N_val, r, t_val, gamma, 'depolarizing', 'II')
        corrs.append(abs(Ic_bw - Ic_rm_val) / k_val)
    corrs = np.array(corrs)
    mask = corrs > 1e-15
    ax.loglog(ts[mask], corrs[mask], 'o-', color=colors_II[ni], ms=4, lw=1,
             label=f'N={N_val}')

# Reference: 1/t
t_ref = np.logspace(0.5, 1.5, 50)
ax.loglog(t_ref, 0.5 / t_ref, 'k:', lw=1.5, alpha=0.6, label=r'$1/t$')

ax.set_xlabel('$t$')
ax.set_ylabel(r'$|I_c^{\rm bw} - I_c^{\rm RM}|/k$')
ax.set_title('(d) Setup II corrections')
ax.legend(fontsize=9)
print(f"Panel (d) done")

fig.suptitle(r'Figure 2: Exact 2D stat-mech, depolarizing noise, $d=2$, $r=1/4$',
             fontsize=13, y=1.01)
fig.tight_layout()
fig.savefig('fig2_exact.pdf', bbox_inches='tight', dpi=150)
fig.savefig('fig2_exact.png', bbox_inches='tight', dpi=150)
print("Saved fig2_exact.pdf/png")
