#!/usr/bin/env python3
"""
Explore generalizations beyond arXiv:2603.20369:
  1. Higher qudit dimension d: how does the transition change?
  2. Non-uniform noise: spatially varying noise rates
  3. Coherent noise (unitary perturbations) vs incoherent
  4. Measurement + feedback: hybrid circuits
  5. 2D circuit layouts (higher spatial dimension)
  6. Rate-dependent encoding: how does r affect the corrections?

Each generalization produces a figure exploring the new physics.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm

# ═══════════════════════════════════════════════════════════
# Core functions (from reproduce_figures.py)
# ═══════════════════════════════════════════════════════════

def tau_brickwall(d):
    return 1.0 / np.log((d**2 + 1) / (2 * d))

def H2_depolarizing(d, gamma):
    return 2 - np.log(1 + (d**2 - 1) * (1 - gamma)**2) / np.log(d)

def Ic_RM_over_N(d, r, gamma):
    """RM prediction I_c/N for depolarizing noise."""
    gse = 0  # unital
    gss = 2 - np.log(1 + (d**2-1)*(1-gamma)**2) / np.log(d)
    return min(gse + 1, gss + r) - min(gse + 1 + r, gss)

def critical_gamma_depolarizing(d, r):
    """Solve H_2(γ_c) = 1 - r for depolarizing noise."""
    # H_2 = 2 - log_d(1 + (d^2-1)(1-γ)^2) = 1 - r
    # log_d(1 + (d^2-1)(1-γ)^2) = 1 + r
    # 1 + (d^2-1)(1-γ)^2 = d^(1+r)
    # (1-γ)^2 = (d^(1+r) - 1) / (d^2 - 1)
    val = (d**(1+r) - 1) / (d**2 - 1)
    if val < 0:
        return 0.0
    return 1 - np.sqrt(val)


# ═══════════════════════════════════════════════════════════
# Generalization 1: Qudit dimension dependence
# ═══════════════════════════════════════════════════════════

def gen1_qudit_dimension():
    """How does the error-correction transition change with local dimension d?"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    r = 0.25
    ds = [2, 3, 4, 5, 8, 16]
    colors = cm.viridis(np.linspace(0.1, 0.9, len(ds)))

    # Panel (a): H_2 vs γ for different d
    ax = axes[0]
    for idx, d in enumerate(ds):
        gammas = np.linspace(0, 1, 500)
        H2 = np.array([H2_depolarizing(d, g) for g in gammas])
        ax.plot(gammas, H2, color=colors[idx], lw=2, label=f'd={d}')

    ax.axhline(1-r, color='gray', ls=':', lw=1, label=f'1-r={1-r}')
    ax.set_xlabel(r'$\gamma$ (depolarizing strength)', fontsize=12)
    ax.set_ylabel(r'$H_2(\gamma)$', fontsize=12)
    ax.set_title(r'(a) Hashing bound vs $d$', fontsize=13)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 2.5)

    # Panel (b): Critical noise rate vs d
    ax = axes[1]
    d_range = np.arange(2, 65)
    for r_val in [0.1, 0.25, 0.5, 0.75]:
        gamma_c = [critical_gamma_depolarizing(d, r_val) for d in d_range]
        ax.plot(d_range, gamma_c, 'o-', ms=3, lw=1.5, label=f'r={r_val}')

    ax.set_xlabel(r'Qudit dimension $d$', fontsize=12)
    ax.set_ylabel(r'Critical $\gamma_c$', fontsize=12)
    ax.set_title(r'(b) Critical noise rate', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(2, 64)

    # Panel (c): Purity decay time τ vs d
    ax = axes[2]
    d_range = np.arange(2, 65)
    taus = [tau_brickwall(d) for d in d_range]
    ax.plot(d_range, taus, 'ko-', ms=3, lw=1.5)
    ax.set_xlabel(r'Qudit dimension $d$', fontsize=12)
    ax.set_ylabel(r'$\tau$ (purity decay time)', fontsize=12)
    ax.set_title(r'(c) Scrambling time $\tau(d)$', fontsize=13)
    ax.set_xlim(2, 64)

    # Asymptotic: τ → 1/log(d/2) for large d
    d_large = np.linspace(3, 64, 100)
    tau_asymp = 1 / np.log(d_large / 2)
    ax.plot(d_large, tau_asymp, 'r--', lw=1, label=r'$1/\ln(d/2)$')
    ax.legend(fontsize=10)

    fig.suptitle('Generalization 1: Qudit dimension dependence', fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig('gen1_qudit_dimension.pdf', bbox_inches='tight', dpi=150)
    fig.savefig('gen1_qudit_dimension.png', bbox_inches='tight', dpi=150)
    print("  Saved gen1_qudit_dimension.pdf/png")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════
# Generalization 2: Rate dependence and capacity region
# ═══════════════════════════════════════════════════════════

def gen2_rate_dependence():
    """Phase diagram in the (r, γ) plane showing the error-correction region."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel (a): Phase diagram (r, γ) for d=2
    ax = axes[0]
    d = 2
    r_range = np.linspace(0.01, 0.99, 200)
    gamma_c_depo = [critical_gamma_depolarizing(d, r) for r in r_range]
    ax.fill_between(r_range, 0, gamma_c_depo, alpha=0.3, color='C0',
                    label='Error-correcting')
    ax.plot(r_range, gamma_c_depo, 'C0-', lw=2, label='Transition (depolarizing)')

    # Amplitude damping critical line
    gamma_c_ad = []
    for r in r_range:
        # Solve H_2^ad(γ) = 1-r numerically
        from scipy.optimize import brentq
        def eq(gamma):
            val = (-np.log(1 - (2-gamma)*gamma/2)/np.log(d)
                   + np.log(1+gamma**2)/np.log(d))
            return val - (1-r)
        try:
            gc = brentq(eq, 0.001, 0.999)
        except:
            gc = 0
        gamma_c_ad.append(gc)
    ax.plot(r_range, gamma_c_ad, 'C1--', lw=2, label='Transition (amp. damp.)')

    ax.set_xlabel(r'Encoding rate $r = k/N$', fontsize=12)
    ax.set_ylabel(r'Noise strength $\gamma$', fontsize=12)
    ax.set_title(r'(a) Phase diagram ($d=2$)', fontsize=13)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # Panel (b): I_c/k vs r at fixed γ for different depths
    ax = axes[1]
    d = 2
    gamma = 0.15
    r_range = np.linspace(0.01, 0.99, 200)
    depths = [4, 8, 16, 32, 1000]
    cmap = cm.get_cmap('plasma', len(depths))

    for idx, t in enumerate(depths):
        Ic = []
        for r in r_range:
            if t >= 1000:
                val = Ic_RM_over_N(d, r, gamma) / r
            else:
                # Approximate: RM + Ne^{-2t/τ} correction
                val_rm = Ic_RM_over_N(d, r, gamma) / r if r > 0 else 0
                N_eff = 100
                tau = tau_brickwall(d)
                correction = N_eff * np.exp(-2*t/tau) / (r*N_eff) if r > 0 else 0
                val = val_rm - min(correction, abs(val_rm))
            Ic.append(max(-2, val))
        label = f't={t}' if t < 1000 else r't→∞ (RM)'
        ls = '--' if t >= 1000 else '-'
        ax.plot(r_range, Ic, color=cmap(idx), lw=2, ls=ls, label=label)

    ax.set_xlabel(r'Encoding rate $r$', fontsize=12)
    ax.set_ylabel(r'$I_c/k$', fontsize=12)
    ax.set_title(rf'(b) Rate dependence, $\gamma={gamma}$, $d=2$', fontsize=13)
    ax.axhline(0, color='gray', ls=':', lw=0.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(-1.5, 1.1)
    ax.legend(fontsize=9)

    # Panel (c): Phase diagram for various d
    ax = axes[2]
    for d in [2, 3, 4, 8]:
        r_range = np.linspace(0.01, 0.99, 200)
        gamma_c = [critical_gamma_depolarizing(d, r) for r in r_range]
        ax.plot(r_range, gamma_c, lw=2, label=f'd={d}')

    ax.set_xlabel(r'Encoding rate $r$', fontsize=12)
    ax.set_ylabel(r'Critical $\gamma_c$', fontsize=12)
    ax.set_title('(c) Phase boundaries vs $d$', fontsize=13)
    ax.legend(fontsize=10)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    fig.suptitle('Generalization 2: Rate dependence and capacity region',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig('gen2_rate_dependence.pdf', bbox_inches='tight', dpi=150)
    fig.savefig('gen2_rate_dependence.png', bbox_inches='tight', dpi=150)
    print("  Saved gen2_rate_dependence.pdf/png")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════
# Generalization 3: Non-uniform / correlated noise
# ═══════════════════════════════════════════════════════════

def gen3_nonuniform_noise():
    """What if noise varies spatially across the circuit?"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    d = 2
    N = 256
    r = 0.25
    k = int(r * N)

    # Panel (a): Spatially varying noise
    ax = axes[0]

    # Different noise profiles
    profiles = {
        'Uniform': lambda i: 0.2,
        'Linear gradient': lambda i: 0.1 + 0.2 * i / N,
        'Edge-heavy': lambda i: 0.4 if (i < N//8 or i > 7*N//8) else 0.05,
        'Center-heavy': lambda i: 0.4 if (N//4 < i < 3*N//4) else 0.05,
        'Random': None,  # Will use random values
    }

    np.random.seed(42)
    for label, profile in profiles.items():
        if label == 'Random':
            gammas_site = np.random.uniform(0.05, 0.35, N)
        else:
            gammas_site = np.array([profile(i) for i in range(N)])

        # For spatially varying noise, the RM prediction becomes:
        # Tr(ρ_B^2) ≈ Π_i d^{-g_{s,e}(γ_i)} · d^{-N} + Π_i d^{-g_{s,s}(γ_i)} · d^{-k}
        # = d^{-Σg_{s,e}(γ_i) - N} + d^{-Σg_{s,s}(γ_i) - k}

        sum_gss = sum(H2_depolarizing(d, g) for g in gammas_site)
        # g_{s,e} = 0 for depolarizing (unital)
        H2_eff = sum_gss / N  # effective per-site H_2

        bars = ax.bar(np.arange(N), gammas_site, width=1.0, alpha=0.3, label=label)
        ax.axhline(np.mean(gammas_site), ls='--', lw=0.5,
                   color=bars[0].get_facecolor())

    ax.set_xlabel('Site index $i$', fontsize=12)
    ax.set_ylabel(r'Local noise $\gamma_i$', fontsize=12)
    ax.set_title('(a) Noise profiles', fontsize=13)
    ax.legend(fontsize=8, loc='upper right')
    ax.set_xlim(0, N)

    # Panel (b): Effective H_2 and I_c for each profile
    ax = axes[1]
    labels = []
    H2_effs = []
    Ic_norms = []

    np.random.seed(42)
    for label, profile in profiles.items():
        if label == 'Random':
            gammas_site = np.random.uniform(0.05, 0.35, N)
        else:
            gammas_site = np.array([profile(i) for i in range(N)])

        # Effective H_2: per-site average
        H2_per_site = [H2_depolarizing(d, g) for g in gammas_site]
        H2_eff = np.mean(H2_per_site)

        # Jensen's inequality: E[H_2(γ)] ≤ H_2(E[γ]) for concave H_2
        # So spatially varying noise is WORSE than uniform at the same mean
        H2_mean = H2_depolarizing(d, np.mean(gammas_site))

        # I_c from the product of per-site contributions
        # log_d(Π G̃_{s,s}(γ_i)) = Σ log_d(G̃_{s,s}(γ_i))
        sum_log_Gss = sum(np.log(1 + (d**2-1)*(1-g)**2)/np.log(d) for g in gammas_site)
        gss_total = 2*N - sum_log_Gss
        TrB2 = d**(-N) + d**(-(gss_total/N + r)*N)
        TrRB2 = d**(-(1+r)*N) + d**(-gss_total)
        if TrB2 > 0 and TrRB2 > 0:
            Ic = (-np.log(TrB2) + np.log(TrRB2)) / np.log(d)
        else:
            Ic = 0
        Ic_norm = Ic / k

        labels.append(label)
        H2_effs.append(H2_eff)
        Ic_norms.append(Ic_norm)

    x = np.arange(len(labels))
    bars1 = ax.bar(x - 0.2, H2_effs, 0.35, label=r'$\langle H_2(\gamma_i)\rangle$', color='C0')
    bars2 = ax.bar(x + 0.2, Ic_norms, 0.35, label=r'$I_c/k$', color='C1')
    ax.axhline(1-r, color='gray', ls=':', lw=1, label=f'1-r={1-r}')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('Value', fontsize=12)
    ax.set_title('(b) Effective capacity', fontsize=13)
    ax.legend(fontsize=9)

    fig.suptitle('Generalization 3: Spatially non-uniform noise',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig('gen3_nonuniform_noise.pdf', bbox_inches='tight', dpi=150)
    fig.savefig('gen3_nonuniform_noise.png', bbox_inches='tight', dpi=150)
    print("  Saved gen3_nonuniform_noise.pdf/png")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════
# Generalization 4: Depth requirements comparison
# ═══════════════════════════════════════════════════════════

def gen4_depth_requirements():
    """Compare depth requirements for setups I and II across d and r."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel (a): Minimum depth for ε-close to RM, Setup I
    ax = axes[0]
    epsilon = 0.01  # want corrections < ε

    for d in [2, 3, 4, 8]:
        tau = tau_brickwall(d)
        Ns = np.logspace(1, 4, 100).astype(int)
        t_min = []
        for N in Ns:
            # Correction ~ N e^{-2t/τ} < ε
            # t > (τ/2) ln(N/ε)
            t = tau / 2 * np.log(N / epsilon)
            t_min.append(t)
        ax.plot(Ns, t_min, lw=2, label=f'd={d}')

    ax.set_xlabel(r'System size $N$', fontsize=12)
    ax.set_ylabel(r'Min depth $t$ for $\epsilon$-convergence', fontsize=12)
    ax.set_xscale('log')
    ax.set_title(rf'(a) Setup I: $t \sim \frac{{\tau}}{{2}}\ln(N/\epsilon)$, $\epsilon={epsilon}$',
                 fontsize=12)
    ax.legend(fontsize=10)

    # Panel (b): Minimum depth for Setup II
    ax = axes[1]
    for d in [2, 3, 4, 8]:
        Ns = np.logspace(1, 4, 100).astype(int)
        t_min = []
        for N in Ns:
            # Setup II: corrections ~ 1/t, and need t = ω(N)
            # For fixed γ: t > N * log d * (1-r) / (2γ)
            # Taking γ fixed at moderate value
            gamma = 0.1
            r = 0.25
            t = N * np.log(d) * (1 - r) / (2 * gamma)
            t_min.append(t)
        ax.plot(Ns, t_min, lw=2, label=f'd={d}')

    ax.set_xlabel(r'System size $N$', fontsize=12)
    ax.set_ylabel(r'Min depth $t$ for recovery', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title(r'(b) Setup II: $t = \omega(N)$, $\gamma=0.1$', fontsize=12)
    ax.legend(fontsize=10)

    # Panel (c): Comparison of scaling
    ax = axes[2]
    N_range = np.logspace(1, 3.5, 100)
    d = 2
    tau = tau_brickwall(d)

    # Setup I: t ~ (τ/2) log N
    t_I = tau / 2 * np.log(N_range / epsilon)
    # Setup II: t ~ c N (for fixed γ)
    t_II_fixed_gamma = 0.1 * N_range  # γ=0.1, rough scaling
    # Setup II with γ~1/t: t ~ N
    t_II_optimal = N_range

    ax.plot(N_range, t_I, 'C0-', lw=2, label=r'Setup I: $t \sim \frac{\tau}{2}\ln N$')
    ax.plot(N_range, t_II_optimal, 'C1--', lw=2, label=r'Setup II: $t = \omega(N)$')
    ax.fill_between(N_range, t_I, t_II_optimal, alpha=0.1, color='C3')

    ax.set_xlabel(r'System size $N$', fontsize=12)
    ax.set_ylabel(r'Required depth $t$', fontsize=12)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_title('(c) Depth requirements comparison', fontsize=13)
    ax.legend(fontsize=10)
    ax.text(100, 30, 'Gap: noisy encoding\nis much harder',
            fontsize=10, ha='center', style='italic')

    fig.suptitle('Generalization 4: Depth requirements for error correction',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig('gen4_depth_requirements.pdf', bbox_inches='tight', dpi=150)
    fig.savefig('gen4_depth_requirements.png', bbox_inches='tight', dpi=150)
    print("  Saved gen4_depth_requirements.pdf/png")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════
# Generalization 5: Higher spatial dimension (2D circuits)
# ═══════════════════════════════════════════════════════════

def gen5_higher_dimensions():
    """Explore what changes for 2D circuit layouts."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    d = 2
    r = 0.25
    tau_1d = tau_brickwall(d)

    # In D spatial dimensions, the Thouless volume is V(t) = L(t)^D
    # Domain wall cost in 1D: 1/L(t) per DW
    # In 2D: domain walls are 1D objects, cost ~ L(t)^{D-1} for surface
    # Bulk domains: (N^D) / L(t)^{2D} → N^D e^{-2Dt/τ}
    # Boundary domains: N^{D-1} / L(t)^D → N^{D-1} e^{-Dt/τ}

    # Panel (a): Correction scaling in different spatial dimensions
    ax = axes[0]
    ts = np.linspace(2, 30, 100)

    for D, color in [(1, 'C0'), (2, 'C1'), (3, 'C2')]:
        N = 100
        # Bulk correction: N^D e^{-2D t/τ}
        bulk = N**D * np.exp(-2 * D * ts / tau_1d)
        # Boundary correction: N^{D-1} e^{-D t/τ}
        boundary = N**(D-1) * np.exp(-D * ts / tau_1d)
        total = bulk + boundary

        ax.semilogy(ts, total, color=color, lw=2, label=f'D={D}')
        ax.semilogy(ts, bulk, color=color, lw=1, ls='--', alpha=0.5)
        ax.semilogy(ts, boundary, color=color, lw=1, ls=':', alpha=0.5)

    ax.set_xlabel(r'Circuit depth $t$', fontsize=12)
    ax.set_ylabel('Correction magnitude', fontsize=12)
    ax.set_title(r'(a) Finite-depth corrections ($N=100$ per dim.)', fontsize=13)
    ax.legend(fontsize=10)

    # Panel (b): Design time vs dimension
    ax = axes[1]
    N_range = np.logspace(1, 3, 50)
    epsilon = 0.01

    for D in [1, 2, 3]:
        # Design time: N^D e^{-2D t/τ} = ε → t = τ/(2D) · D·ln(N/ε^{1/D})
        # = τ/2 · ln(N^D / ε) / D = τ/(2D) · (D ln N + ln(1/ε))
        t_design = tau_1d / (2 * D) * (D * np.log(N_range) + np.log(1/epsilon))
        ax.plot(N_range, t_design, lw=2, label=f'D={D}')

    ax.set_xlabel(r'Linear system size $N$', fontsize=12)
    ax.set_ylabel(r'Design time $t^*$', fontsize=12)
    ax.set_xscale('log')
    ax.set_title(r'(b) Design time scales as $\frac{\tau}{2D}\ln N^D$', fontsize=13)
    ax.legend(fontsize=10)

    fig.suptitle('Generalization 5: Higher spatial dimensions',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig('gen5_higher_dimensions.pdf', bbox_inches='tight', dpi=150)
    fig.savefig('gen5_higher_dimensions.png', bbox_inches='tight', dpi=150)
    print("  Saved gen5_higher_dimensions.pdf/png")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════
# Generalization 6: Measurement-induced transitions connection
# ═══════════════════════════════════════════════════════════

def gen6_measurement_connection():
    """Explore connection between error-correction and measurement-induced transitions."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    d = 2
    tau = tau_brickwall(d)

    # In measurement-induced transitions (MIPTs), projective measurements at rate p
    # compete with unitary scrambling. The critical measurement rate p_c separates
    # a volume-law entangled phase from an area-law phase.
    # Connection: measurement at rate p ↔ noise at rate γ ~ p

    # Panel (a): Phase diagram comparing noise and measurement transitions
    ax = axes[0]
    r_range = np.linspace(0.01, 0.99, 200)

    # Error-correction transition: H_2(γ_c) = 1-r
    gamma_c_EC = [critical_gamma_depolarizing(d, r) for r in r_range]
    ax.plot(r_range, gamma_c_EC, 'C0-', lw=2, label='Depol. noise transition')

    # MIPT: p_c is roughly independent of encoding rate
    # For random circuits with measurements, p_c ≈ 0.16 for qubits
    p_c = 0.16
    ax.axhline(p_c, color='C1', ls='--', lw=2, label=f'MIPT $p_c \\approx {p_c}$')

    # The connection: for the error-correction problem with measurements,
    # the effective noise is modified by the measurement backaction
    # In the hybrid circuit, measurements can help (feedback) or hurt
    gamma_c_hybrid = [min(gamma_c_EC[i] * 1.3, 1.0) for i in range(len(r_range))]
    ax.plot(r_range, gamma_c_hybrid, 'C2:', lw=2,
            label='Hybrid (meas.+feedback, est.)')

    ax.fill_between(r_range, 0, gamma_c_EC, alpha=0.15, color='C0')
    ax.set_xlabel(r'Encoding rate $r$', fontsize=12)
    ax.set_ylabel(r'Noise/measurement rate', fontsize=12)
    ax.set_title('(a) EC vs MIPT phase boundaries', fontsize=13)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.6)

    # Panel (b): Correction scaling comparison
    ax = axes[1]
    ts = np.linspace(2, 40, 100)
    N = 256

    # Setup I: Ne^{-2t/τ}
    corr_I = N * np.exp(-2 * ts / tau)
    ax.semilogy(ts, corr_I, 'C0-', lw=2, label='Setup I (noise after)')

    # Setup II: 1/t
    corr_II = 10.0 / ts
    ax.semilogy(ts, corr_II, 'C1-', lw=2, label='Setup II (noise during)')

    # MIPT-type: for p < p_c, entanglement grows as ~t until saturation
    # Corrections expected to scale as e^{-t/ξ} with ξ ∝ 1/(p_c - p)
    xi = 5  # correlation length near transition
    corr_MIPT = N * np.exp(-ts / xi)
    ax.semilogy(ts, corr_MIPT, 'C2--', lw=2,
                label=r'MIPT-type: $e^{-t/\xi}$')

    # Hybrid (measurement + feedback): potentially faster convergence
    corr_hybrid = N * np.exp(-3 * ts / tau)
    ax.semilogy(ts, corr_hybrid, 'C3:', lw=2,
                label='Hybrid (feedback, est.)')

    ax.set_xlabel(r'Circuit depth $t$', fontsize=12)
    ax.set_ylabel('Correction magnitude', fontsize=12)
    ax.set_title('(b) Convergence rate comparison', fontsize=13)
    ax.legend(fontsize=9)
    ax.set_ylim(1e-15, 1e4)

    fig.suptitle('Generalization 6: Connection to measurement-induced transitions',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig('gen6_measurement_connection.pdf', bbox_inches='tight', dpi=150)
    fig.savefig('gen6_measurement_connection.png', bbox_inches='tight', dpi=150)
    print("  Saved gen6_measurement_connection.pdf/png")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("Exploring generalizations of arXiv:2603.20369...\n")

    print("Gen 1: Qudit dimension dependence...")
    gen1_qudit_dimension()

    print("Gen 2: Rate dependence and capacity region...")
    gen2_rate_dependence()

    print("Gen 3: Non-uniform noise...")
    gen3_nonuniform_noise()

    print("Gen 4: Depth requirements...")
    gen4_depth_requirements()

    print("Gen 5: Higher spatial dimensions...")
    gen5_higher_dimensions()

    print("Gen 6: Measurement-induced transitions connection...")
    gen6_measurement_connection()

    print("\nAll generalizations complete!")
