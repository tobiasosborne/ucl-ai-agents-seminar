# 05 — Live Demo: FeO Potential Energy Curves

An audience member at the seminar challenged the AI coding agent to
tackle a real computational-physics problem on the spot:

> **"Compute the potential energy curves for the lowest ~25 electronic
> states of FeO (iron monoxide)."**

FeO is a notoriously difficult molecule for electronic-structure theory.
Its open 3d shell, strong multi-reference character, and dense manifold of
near-degenerate quintet/triplet/septet states make it a genuine test of
both scientific knowledge and software-engineering ability.

## What the agent produced (in ~20 minutes)

Starting from an empty directory, the agent:

1. Created a Python virtual environment and installed PySCF.
2. Wrote and ran a **smoke test** (`smoke_test_feo.py`) to validate the
   full ROHF -> SA-CASSCF -> CASCI -> SC-NEVPT2 pipeline at a single
   geometry.
3. Wrote a **production script** (`compute_feo_pec.py`) that scans ~35
   bond lengths and computes ~25 electronic states at each geometry using
   SA-CASSCF(12,12) + SC-NEVPT2 with X2C relativistic corrections and
   Douglas-Kroll basis sets.
4. Wrote an **analysis script** (`plot_feo_pec.py`) to plot the PECs and
   extract spectroscopic constants (r_e, T_e, omega_e, D_e).
5. Diagnosed and fixed a PySCF-specific bug: NEVPT2 cannot be called
   directly on a state-averaged CASSCF object — an intermediate multi-root
   CASCI step is required (cf. PySCF `examples/mrpt/41-for_state_average.py`).

The smoke test completed successfully. The full production run is too
expensive for a laptop (~hours per geometry at the (12,12) active space)
but the methodology is validated end-to-end.

## Files

| File | Purpose |
|------|---------|
| `smoke_test_feo.py` | Minimal single-point test: 1 geometry, small active space (8e,8o), 2 roots |
| `compute_feo_pec.py` | Production PEC scan with SA-CASSCF(12,12) + SC-NEVPT2 |
| `plot_feo_pec.py` | Plot curves, identify degenerate pairs, extract spectroscopic constants |
| `hello.py` | Initial venv sanity check |

## Method summary

- **Active space:** 12 electrons in 12 orbitals (Fe 3d + 4s + O 2p + correlating 4d-like)
- **Dynamic correlation:** Strongly-contracted NEVPT2
- **Relativistic treatment:** Exact two-component (X2C) scalar Hamiltonian
- **Basis sets:** cc-pVTZ-DK (production) / cc-pVDZ-DK (quick test)
- **States:** Quintets (13 roots), triplets (11 roots), septets (7 roots) across C2v irreps

## Running

```bash
# Create venv and install dependencies
python3 -m venv feo_venv
feo_venv/bin/pip install pyscf numpy scipy matplotlib h5py

# Quick validation (single geometry, ~3 min)
feo_venv/bin/python3 smoke_test_feo.py

# Reduced test (5 geometries, quintets only, ~2 hours)
feo_venv/bin/python3 compute_feo_pec.py --quick --reduced

# Full production (~days on a workstation, suited for HPC)
feo_venv/bin/python3 compute_feo_pec.py

# Plot results
feo_venv/bin/python3 plot_feo_pec.py
```

## Point of the demo

This exercise demonstrated that an AI coding agent can:

- Identify the correct multi-reference quantum chemistry methodology
  for a genuinely hard open-shell transition-metal problem
- Navigate PySCF's API and work around non-obvious library constraints
- Produce a complete, runnable computational pipeline — not just code
  snippets — in minutes rather than the hours or days it would take
  a researcher starting from scratch
