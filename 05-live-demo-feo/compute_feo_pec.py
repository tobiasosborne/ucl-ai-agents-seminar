#!/usr/bin/env python3
"""
Potential Energy Curves of FeO — Lowest ~25 Electronic States

Method:  SA-CASSCF(12,12) + SC-NEVPT2
Basis:   cc-pVTZ-DK
Relativ: Exact 2-Component (X2C) scalar relativistic Hamiltonian

Active space (12 electrons in 12 orbitals):
  Fe 3d (5) + Fe 4s (1) + O 2p (3) + correlating 4d-like (3)

States computed in C2v subgroup of C-inf-v:
  Quintets  (2S+1=5): ground-state manifold, X ^5Delta, etc.
  Triplets  (2S+1=3): low-lying excited states
  Septets   (2S+1=7): high-spin states from Fe(^5D/^7S) + O(^3P)

Usage:
  python compute_feo_pec.py              # full calculation
  python compute_feo_pec.py --quick      # test run (5 geometries, cc-pVDZ-DK)
  python compute_feo_pec.py --resume     # resume interrupted calculation
"""

import numpy as np
import json
import os
import sys
import time
import argparse
from pathlib import Path


# ===================================================================
# Configuration
# ===================================================================

NCAS = 12       # active orbitals
NELECAS = 12    # active electrons

# States: (multiplicity, C2v_irrep, nroots)
# Enough to capture ~25 unique Λ-S states after identifying degenerate pairs
STATES_FULL = [
    # Quintets — dominant low-energy manifold
    (5, 'A1', 4),   # Sigma+, Delta(a), Gamma(a), ...
    (5, 'A2', 3),   # Sigma-, Delta(b), Gamma(b)
    (5, 'B1', 3),   # Pi(a), Phi(a), ...
    (5, 'B2', 3),   # Pi(b), Phi(b), ...
    # Triplets
    (3, 'A1', 3),
    (3, 'A2', 2),
    (3, 'B1', 3),
    (3, 'B2', 3),
    # Septets
    (7, 'A1', 2),
    (7, 'A2', 1),
    (7, 'B1', 2),
    (7, 'B2', 2),
]

# Reduced set: quintets only, fewer roots — ~10 states, much faster
STATES_REDUCED = [
    (5, 'A1', 3),   # Sigma+, Delta(a), ...
    (5, 'A2', 2),   # Sigma-, Delta(b)
    (5, 'B1', 2),   # Pi(a), Phi(a)
    (5, 'B2', 2),   # Pi(b), Phi(b)
]


def make_r_grid(quick=False):
    """Bond length grid in Angstrom."""
    if quick:
        return np.array([1.45, 1.62, 1.80, 2.20, 3.00])
    r = np.concatenate([
        np.arange(1.30, 1.55, 0.05),
        np.arange(1.55, 1.75, 0.02),
        np.arange(1.75, 2.10, 0.05),
        np.arange(2.10, 3.60, 0.10),
    ])
    return np.sort(np.unique(np.round(r, 4)))


def build_mol(r_ang, spin_2S, basis):
    """Build FeO Mole object at bond length r (Angstrom) with given 2S."""
    from pyscf import gto
    mol = gto.Mole()
    mol.atom = f'Fe 0 0 0; O 0 0 {r_ang}'
    mol.basis = basis
    mol.spin = spin_2S
    mol.charge = 0
    mol.symmetry = 'C2v'
    mol.unit = 'Angstrom'
    mol.max_memory = 8000
    mol.verbose = 3
    mol.build()
    return mol


def run_rohf(mol):
    """ROHF with X2C scalar relativistic Hamiltonian."""
    from pyscf import scf
    mf = scf.ROHF(mol).x2c()
    mf.max_cycle = 300
    mf.level_shift = 0.3
    mf.conv_tol = 1e-9
    mf.kernel()
    if not mf.converged:
        # Retry with ATOM initial guess and stronger shift
        mf.level_shift = 0.8
        mf.init_guess = 'atom'
        mf.kernel()
    return mf


def run_casscf_nevpt2(mf, ncas, nelecas, irrep, nroots, prev_mo=None):
    """
    State-averaged CASSCF + SC-NEVPT2 for a given C2v irrep.

    Returns dict with energies, convergence flag, and MO coefficients.
    """
    from pyscf import mcscf
    from pyscf.mrpt import nevpt2 as nevpt2_mod

    mc = mcscf.CASSCF(mf, ncas, nelecas)
    mc.fcisolver.wfnsym = irrep
    mc.fcisolver.max_cycle = 200
    mc.max_cycle_macro = 200
    mc.max_cycle_micro = 10
    mc.conv_tol = 1e-8
    mc.natorb = True

    if nroots > 1:
        mc = mcscf.state_average_(mc, [1.0 / nroots] * nroots)

    # Use previous MOs for orbital continuity along PEC
    if prev_mo is not None:
        mc.kernel(prev_mo)
    else:
        mc.kernel()

    # Extract CASSCF energies
    if nroots > 1:
        casscf_energies = list(mc.e_states)
    else:
        casscf_energies = [mc.e_tot]

    # For NEVPT2 on SA-CASSCF: PySCF requires a multi-root CASCI step
    # using the state-averaged orbitals.
    mc_ci = mcscf.CASCI(mf, ncas, nelecas)
    mc_ci.fcisolver.wfnsym = irrep
    mc_ci.fcisolver.nroots = nroots
    mc_ci.kernel(mc.mo_coeff)

    nevpt2_energies = []
    for iroot in range(nroots):
        try:
            e_corr = nevpt2_mod.NEVPT(mc_ci, root=iroot).kernel()
            e_casci = mc_ci.e_tot[iroot] if nroots > 1 else mc_ci.e_tot
            nevpt2_energies.append(e_casci + e_corr)
        except Exception as e:
            print(f"    NEVPT2 root {iroot} failed: {e}")
            nevpt2_energies.append(None)

    return {
        'casscf': casscf_energies,
        'nevpt2': nevpt2_energies,
        'converged': bool(mc.converged),
        'mo_coeff': mc.mo_coeff,
    }


def main():
    parser = argparse.ArgumentParser(description='FeO PEC computation')
    parser.add_argument('--quick', action='store_true',
                        help='Quick test: 5 geometries, cc-pVDZ-DK basis')
    parser.add_argument('--resume', action='store_true',
                        help='Resume from saved intermediate results')
    parser.add_argument('--reduced', action='store_true',
                        help='Reduced scope: quintets only, fewer roots (~10 states)')
    args = parser.parse_args()

    basis = 'cc-pvdz-dk' if args.quick else 'cc-pvtz-dk'
    r_grid = make_r_grid(quick=args.quick)
    states = STATES_REDUCED if args.reduced else STATES_FULL

    outdir = Path('feo_pec_results')
    outdir.mkdir(exist_ok=True)
    result_file = outdir / 'pec_data.json'

    # Load previous results if resuming
    all_results = {}
    if args.resume and result_file.exists():
        with open(result_file) as f:
            saved = json.load(f)
            all_results = saved.get('results', {})
        print(f"Resuming: loaded {len(all_results)} bond lengths from previous run")

    # Track MOs for orbital continuity (keyed by (mult, irrep))
    prev_mos = {}

    total_geoms = len(r_grid)
    total_blocks = len(states)
    t_start = time.time()

    for ir, r in enumerate(r_grid):
        r_key = f"{r:.4f}"

        # Check if this geometry is already complete
        if r_key in all_results:
            existing = all_results[r_key]
            if all(f"{m}_{sym}" in existing for m, sym, _ in states):
                print(f"[{ir+1}/{total_geoms}] r={r:.4f} A -- already done, skipping")
                continue

        print(f"\n{'='*65}")
        print(f"  r = {r:.4f} A   [{ir+1}/{total_geoms}]")
        print(f"{'='*65}")

        if r_key not in all_results:
            all_results[r_key] = {}

        # Group states by spin to reuse ROHF
        from collections import defaultdict
        spin_groups = defaultdict(list)
        for mult, irrep, nroots in states:
            spin_groups[mult].append((irrep, nroots))

        for mult in sorted(spin_groups.keys()):
            spin_2S = mult - 1

            # Check if all irreps for this spin already computed
            irreps_todo = []
            for irrep, nroots in spin_groups[mult]:
                state_key = f"{mult}_{irrep}"
                if state_key in all_results[r_key]:
                    continue
                irreps_todo.append((irrep, nroots))

            if not irreps_todo:
                continue

            # Single ROHF per spin multiplicity
            print(f"\n  ROHF: 2S+1 = {mult}")
            try:
                mol = build_mol(r, spin_2S, basis)
                mf = run_rohf(mol)
                if not mf.converged:
                    print(f"  WARNING: ROHF not converged for 2S+1={mult}")
            except Exception as e:
                print(f"  ROHF FAILED for 2S+1={mult}: {e}")
                for irrep, nroots in irreps_todo:
                    state_key = f"{mult}_{irrep}"
                    all_results[r_key][state_key] = {
                        'casscf': [None] * nroots,
                        'nevpt2': [None] * nroots,
                        'converged': False,
                        'error': f'ROHF failed: {e}',
                    }
                continue

            # CASSCF + NEVPT2 for each irrep
            for irrep, nroots in irreps_todo:
                state_key = f"{mult}_{irrep}"
                print(f"\n  --- CASSCF({NELECAS},{NCAS}) "
                      f"2S+1={mult} {irrep} nroots={nroots} ---")
                t0 = time.time()

                try:
                    result = run_casscf_nevpt2(
                        mf, NCAS, NELECAS, irrep, nroots,
                        prev_mo=prev_mos.get((mult, irrep))
                    )
                    prev_mos[(mult, irrep)] = result['mo_coeff']

                    all_results[r_key][state_key] = {
                        'casscf': [float(e) for e in result['casscf']],
                        'nevpt2': [float(e) if e is not None else None
                                   for e in result['nevpt2']],
                        'converged': result['converged'],
                    }
                    dt = time.time() - t0
                    conv = "OK" if result['converged'] else "NOT CONVERGED"
                    print(f"  {conv} | {dt:.1f}s | "
                          f"E(CAS)={result['casscf'][0]:.6f}")

                except Exception as e:
                    dt = time.time() - t0
                    print(f"  FAILED ({dt:.1f}s): {e}")
                    all_results[r_key][state_key] = {
                        'casscf': [None] * nroots,
                        'nevpt2': [None] * nroots,
                        'converged': False,
                        'error': str(e),
                    }

        # Save after each geometry
        with open(result_file, 'w') as f:
            json.dump({
                'config': {
                    'method': 'SA-CASSCF + SC-NEVPT2',
                    'basis': basis,
                    'ncas': NCAS,
                    'nelecas': NELECAS,
                    'relativistic': 'X2C',
                    'states': [[m, s, n] for m, s, n in states],
                    'r_grid': r_grid.tolist(),
                },
                'results': all_results,
            }, f, indent=2)

        elapsed = time.time() - t_start
        print(f"\n  Saved. Elapsed: {elapsed/60:.1f} min")

    total_time = time.time() - t_start
    print(f"\n\nDone! Total time: {total_time/3600:.2f} hours")
    print(f"Results: {result_file}")


if __name__ == '__main__':
    main()
