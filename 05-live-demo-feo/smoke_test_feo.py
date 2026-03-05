#!/usr/bin/env python3
"""Minimal smoke test: 1 geometry, 1 spin, 1 irrep, small active space."""
import time
t0 = time.time()

from pyscf import gto, scf, mcscf
from pyscf.mrpt import nevpt2

mol = gto.Mole()
mol.atom = 'Fe 0 0 0; O 0 0 1.62'
mol.basis = 'cc-pvdz-dk'
mol.spin = 4   # quintet: 2S=4
mol.symmetry = 'C2v'
mol.max_memory = 4000
mol.verbose = 4
mol.build()

print(f"\nAOs: {mol.nao_nr()}, Symmetry: {mol.groupname}")

# ROHF + X2C
mf = scf.ROHF(mol).x2c()
mf.max_cycle = 300
mf.level_shift = 0.3
mf.kernel()
print(f"\nROHF converged: {mf.converged}, E = {mf.e_tot:.8f}")

# SA-CASSCF (8e, 8o)
ncas, nelecas, nroots = 8, 8, 2
mc = mcscf.CASSCF(mf, ncas, nelecas)
mc.fcisolver.wfnsym = 'A1'
mc.max_cycle_macro = 100
mc.conv_tol = 1e-7
mc = mcscf.state_average_(mc, [1.0/nroots]*nroots)
mc.kernel()

print(f"\nCASSCF converged: {mc.converged}")
print(f"SA-CASSCF energies: {mc.e_states}")

# For NEVPT2 on SA-CASSCF: need multi-root CASCI with the SA orbitals
mc_nevpt = mcscf.CASCI(mf, ncas, nelecas)
mc_nevpt.fcisolver.wfnsym = 'A1'
mc_nevpt.fcisolver.nroots = nroots
mc_nevpt.kernel(mc.mo_coeff)

print(f"\nCASCI energies (with SA orbitals): {mc_nevpt.e_tot}")

for iroot in range(nroots):
    e_corr = nevpt2.NEVPT(mc_nevpt, root=iroot).kernel()
    e_total = mc_nevpt.e_tot[iroot] + e_corr
    print(f"NEVPT2 root {iroot}: E_corr = {e_corr:.8f}, E_total = {e_total:.8f}")

dt = time.time() - t0
print(f"\nSmoke test completed in {dt:.1f} s")
