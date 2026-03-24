#!/usr/bin/env python3
"""
Exact 2D stat-mech computation for brickwall Haar circuits.

KEY INSIGHT: The permutation states |e⟩⟩, |s⟩⟩ are NON-ORTHOGONAL
(⟨⟨e|s⟩⟩ = d). The gate tensor must be expressed in the COEFFICIENT
basis: C = G_2site^{-1} × T_overlap, where G_2site is the 2-site Gram
matrix and T_overlap is the raw overlap transfer matrix.

The partition function is:
  Z = Σ_σ c_final(σ) × Π_j top_overlap_j(σ_j)

where c is evolved by C at each gate, and top_overlap is the physical
overlap of the top boundary with the basis states.
"""
import numpy as np


# ═══════════════════════════════════════════════════════════
# Fundamental objects
# ═══════════════════════════════════════════════════════════

def gram_d(d):
    """Single-site Gram matrix G[π,σ] = ⟨⟨π|σ⟩⟩ for {e,s} on d-dim."""
    return np.array([[d**2, d], [d, d**2]], dtype=np.float64)


def noisy_gram_depolarizing(d, gamma):
    """G̃[π,σ] = ⟨⟨π|N^{⊗2}|σ⟩⟩ for depolarizing noise."""
    return np.array([
        [d**2, d],
        [d, 1 + (d**2 - 1) * (1 - gamma)**2]
    ], dtype=np.float64)


def gate_coefficient_matrix(d):
    """4×4 coefficient transfer matrix C for a 2-site Haar gate.

    C = G_2site^{-1} × T_overlap

    where T_overlap is the raw overlap transfer matrix and
    G_2site = G_d ⊗ G_d is the 2-site Gram matrix.

    Basis order: (e,e)=0, (e,s)=1, (s,e)=2, (s,s)=3
    """
    q = d**2
    # Weingarten on q-dim space
    Wg = np.linalg.inv(np.array([[q**2, q], [q, q**2]], dtype=np.float64))
    # Single-site overlaps
    G = gram_d(d)

    # Raw overlap transfer matrix T[(σ'₁,σ'₂),(σ₁,σ₂)]
    configs = [(0, 0), (0, 1), (1, 0), (1, 1)]
    T_raw = np.zeros((4, 4))
    for oi, (s1o, s2o) in enumerate(configs):
        for ii, (s1i, s2i) in enumerate(configs):
            for pi in range(2):
                for rho in range(2):
                    T_raw[oi, ii] += (
                        Wg[pi, rho]
                        * G[s1o, pi] * G[s2o, pi]
                        * G[rho, s1i] * G[rho, s2i]
                    )

    # 2-site Gram matrix
    G2 = np.kron(G, G)
    # Coefficient matrix
    C = np.linalg.solve(G2, T_raw)
    return C


def ancilla_coefficients(d):
    """Expansion of |0,0⟩^{⊗2} in the {|e⟩⟩, |s⟩⟩} basis.

    |0,0⟩^{⊗2} = c_e |e⟩⟩ + c_s |s⟩⟩ where:
    G @ [c_e, c_s]^T = [1, 1]^T  (since ⟨⟨e|00⟩⟩ = ⟨⟨s|00⟩⟩ = 1)
    """
    G = gram_d(d)
    return np.linalg.solve(G, np.array([1.0, 1.0]))


# ═══════════════════════════════════════════════════════════
# Temporal transfer matrix (exact, N ≤ ~22)
# ═══════════════════════════════════════════════════════════

def _apply_gate_coeff(state, C_gate, N, si, sj):
    """Apply coefficient gate matrix C to sites (si, sj).

    state: array of shape (2,)*N — coefficients in {e,s} basis
    C_gate: 4×4 coefficient matrix, basis (ee,es,se,ss)
    """
    C4 = C_gate.reshape(2, 2, 2, 2)  # [out_i, out_j, in_i, in_j]
    ndim = N
    in_labels = list(range(ndim))
    new_i, new_j = ndim, ndim + 1
    out_labels = list(in_labels)
    out_labels[si] = new_i
    out_labels[sj] = new_j
    t4_labels = [new_i, new_j, in_labels[si], in_labels[sj]]
    return np.einsum(C4, t4_labels, state, in_labels, out_labels, optimize=True)


def purity_exact(d, N, k, t, gamma=0.0, noise_type='depolarizing',
                 setup='I', which='B'):
    """Compute E_U[Tr(ρ_{B or RB}^2)] exactly via temporal transfer matrix.

    Feasible for N ≤ ~22.
    """
    C_gate = gate_coefficient_matrix(d)
    G = gram_d(d)
    c_anc = ancilla_coefficients(d)  # [c_e, c_s] for ancilla sites

    if noise_type == 'depolarizing':
        G_noisy = noisy_gram_depolarizing(d, gamma)
    else:
        raise NotImplementedError(f"noise_type={noise_type}")

    # --- Build noise-modified gate for Setup II ---
    if setup == 'II' and gamma > 0:
        # Noise after each gate modifies the overlap structure.
        # Per-site noise in the replicated picture: N^{⊗2} acts on each site.
        # In coefficient basis: the noise transforms coefficients via
        # c → G^{-1} × G̃ × c  (where G̃ is the noisy Gram matrix).
        # For 2 sites: noise_coeff = (G2^{-1} × G̃_2) applied after gate.
        G_inv = np.linalg.inv(G)
        noise_1site = G_inv @ G_noisy  # coefficient-basis noise per site
        noise_2site = np.kron(noise_1site, noise_1site)  # 4×4
        C_gate_noisy = noise_2site @ C_gate
    else:
        C_gate_noisy = C_gate

    # --- Initialize bottom boundary (coefficients) ---
    state = np.ones((2,) * N, dtype=np.float64)

    for site in range(N):
        if site < k:
            # R-site: |v_B⟩⟩ = |e⟩⟩ for Tr(ρ_B^2), |s⟩⟩ for Tr(ρ_{RB}^2)
            if which == 'B':
                coeff = np.array([1.0, 0.0])  # |e⟩⟩
            else:
                coeff = np.array([0.0, 1.0])  # |s⟩⟩
        else:
            # Ancilla: |0,0⟩^{⊗2}
            coeff = c_anc.copy()

        slc = [slice(None)] * N
        for sigma in range(2):
            slc_s = list(slc)
            slc_s[site] = sigma
            state[tuple(slc_s)] *= coeff[sigma]

    state *= d**(-2 * k)  # normalization from |v_{B/RB}⟩⟩

    # --- Apply brickwall layers ---
    C_use = C_gate_noisy if (setup == 'II' and gamma > 0) else C_gate

    for layer in range(t):
        offset = layer % 2
        for gs in range(offset, N - 1, 2):
            state = _apply_gate_coeff(state, C_use, N, gs, gs + 1)

    # --- Top boundary (overlap) ---
    # Z = Σ_σ c(σ) × Π_j top_j(σ_j)
    # For Setup I: top_j(σ) = G̃[1, σ] = noisy overlap with swap
    # For Setup II: top_j(σ) = G[1, σ] = noiseless overlap with swap
    if setup == 'I':
        top_per_site = G_noisy[1, :]  # [G̃_{s,e}, G̃_{s,s}]
    else:
        top_per_site = G[1, :]  # [G_{s,e}, G_{s,s}] = [d, d²]

    result = state.copy()
    for site in range(N):
        slc = [slice(None)] * N
        for sigma in range(2):
            slc_s = list(slc)
            slc_s[site] = sigma
            result[tuple(slc_s)] *= top_per_site[sigma]

    return np.sum(result)


# ═══════════════════════════════════════════════════════════
# High-level interface
# ═══════════════════════════════════════════════════════════

def coherent_info(d, N, r, t, gamma=0.0, noise_type='depolarizing', setup='I'):
    """Compute annealed 2-Rényi coherent information I_c."""
    k = int(r * N)
    TrB2 = purity_exact(d, N, k, t, gamma, noise_type, setup, 'B')
    TrRB2 = purity_exact(d, N, k, t, gamma, noise_type, setup, 'RB')
    if TrB2 <= 0 or TrRB2 <= 0:
        return float('nan')
    return -np.log(TrB2) / np.log(d) + np.log(TrRB2) / np.log(d)


def Ic_RM(d, N, r, gamma, noise_type='depolarizing'):
    """Random-matrix (t→∞) prediction for I_c."""
    G = gram_d(d)
    Gn = noisy_gram_depolarizing(d, gamma) if noise_type == 'depolarizing' else None
    g = np.log(G) / np.log(d) - np.log(np.maximum(Gn, 1e-300)) / np.log(d)
    gse, gss = g[1, 0], g[1, 1]
    TrB2 = d**(-(gse + 1) * N) + d**(-(gss + r) * N)
    TrRB2 = d**(-(gse + 1 + r) * N) + d**(-gss * N)
    if TrB2 <= 0 or TrRB2 <= 0:
        return float('nan')
    return -np.log(TrB2) / np.log(d) + np.log(TrRB2) / np.log(d)


def H2_depolarizing(d, gamma):
    return 2 - np.log(1 + (d**2 - 1) * (1 - gamma)**2) / np.log(d)


def tau_brickwall(d):
    return 1.0 / np.log((d**2 + 1) / (2 * d))


# ═══════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════

if __name__ == '__main__':
    d = 2

    print("=== Gate coefficient matrix (d=2) ===")
    C = gate_coefficient_matrix(d)
    print(C.round(6))
    print(f"Rank: {np.linalg.matrix_rank(C, tol=1e-10)}")
    print(f"Eigenvalues: {np.linalg.eigvals(C).real.round(6)}")
    print()

    print("=== Ancilla coefficients ===")
    ca = ancilla_coefficients(d)
    print(f"c_e = c_s = {ca[0]:.6f} (expected {1/(d*(d+1))})")
    print()

    print("=== Validation: t=0 (identity circuit) ===")
    for N in [4, 6, 8]:
        k = N // 4
        TrB = purity_exact(d, N, k, 0, 0.0, 'depolarizing', 'I', 'B')
        TrRB = purity_exact(d, N, k, 0, 0.0, 'depolarizing', 'I', 'RB')
        Ic = -np.log(TrB)/np.log(d) + np.log(TrRB)/np.log(d)
        print(f"  N={N}, k={k}: Tr(ρ_B²)={TrB:.6f} (expect {d**(-k):.6f}), "
              f"Tr(ρ_RB²)={TrRB:.6f} (expect 1), Ic={Ic:.4f} (expect {k})")
    print()

    print("=== Validation: convergence to RM ===")
    for N in [4, 8, 12]:
        k = N // 4
        for gamma in [0.0, 0.2]:
            Ic_rm = Ic_RM(d, N, 0.25, gamma)
            for t in [2, 6, 12, 20]:
                Ic_bw = coherent_info(d, N, 0.25, t, gamma)
                diff = Ic_bw - Ic_rm
                print(f"  N={N:2d} γ={gamma} t={t:2d}: Ic_bw={Ic_bw:8.4f} "
                      f"Ic_RM={Ic_rm:8.4f} diff={diff:+.2e}")
            print()

    print("=== Setup II validation ===")
    for N in [4, 8]:
        k = N // 4
        for gamma in [0.1]:
            for t in [2, 6, 12]:
                Ic_bw = coherent_info(d, N, 0.25, t, gamma, 'depolarizing', 'II')
                print(f"  N={N:2d} γ={gamma} t={t:2d}: Ic_bw={Ic_bw:8.4f}")
        print()
