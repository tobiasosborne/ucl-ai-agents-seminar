# Handoff: Formalization of arXiv:2603.20369

## What was done

Complete formalization, verification, figure reproduction, and generalization study of
"Error-Correction Transitions in Finite-Depth Quantum Channels" (Sauliere, Lami, Ribeiro,
De Luca, De Nardis, 2026).

### Paper download
- PDF and TeX source downloaded from arXiv via `playwright-cli` and `curl`
- Local copy at `paper/source/main.tex` (ground truth for all equation references)

### Formalization (`formalization/`)
- **27 proof nodes**, all validated by adversarial verifier (0 challenges)
- 8 definitions, 1 external reference, append-only ledger (105 events)
- `af handoff` confirms: 100% complete, no open work
- Ground truth: 13 key equations verified by string match against TeX source
- Computational verification: independent Python script confirms all analytical formulas

### Figure reproduction (`figures/`)
- `reproduce_figures.py` — approximate 1D domain-wall model (qualitative, fast)
- `exact_statmech.py` — **exact 2D stat-mech** using coefficient transfer matrix
- `fig2_exact.py` — correct Figure 2 using exact computation (N≤20)
- `generalizations.py` — 6 generalization studies

### Report (`report/`)
- `report.tex` / `report.pdf` — 16-page pdflatex report referencing all 27 proof nodes

## Key technical finding

The permutation states {|e⟩⟩, |s⟩⟩} are **non-orthogonal** (⟨⟨e|s⟩⟩ = d). The gate
transfer matrix must use the **coefficient basis**: `C = G_2site⁻¹ × T_overlap`, not
the raw overlap matrix. This was the root cause of incorrect results in the first
reproduction attempt. Fixed in `exact_statmech.py`.

## What remains (next session)

### Priority 1: Spatial transfer matrix for N=512
The exact temporal TM (2^N states) works for N ≤ 20. The paper's Figure 2 uses N=512.
To match exactly, implement the **spatial transfer matrix**:
- State = temporal config of one column, dimension 2^(t/2)
- For t=30: dim = 2^15 = 32K — very feasible
- Apply N times along the chain → works for any N
- Core idea: contract the 2D partition function column-by-column instead of row-by-row

### Priority 2: Setup II noise-in-gate refinement
The current Setup II implementation (`noise_coeff = G⁻¹ G̃` per site folded into gate)
gives correct qualitative behavior but needs validation against the paper's exact boundary
condition treatment (Eq. 19-20). The noise appears in **both replicas** and modifies
both top and bottom boundary conditions differently from Setup I.

### Priority 3: Update report with exact figures
Replace the approximate figures in `report/report.tex` with the exact ones, and add
discussion of the non-orthogonal basis issue.

## File inventory

```
demoform/
├── paper/
│   ├── paper.pdf                    # Downloaded paper
│   ├── source/main.tex              # TeX source (ground truth)
│   └── source/figures/              # Paper's original figures
├── formalization/
│   ├── meta.json                    # af proof metadata
│   ├── ledger/                      # Append-only proof ledger (105 events)
│   ├── externals/                   # External references
│   ├── proof_tree.tex               # af export --format latex
│   ├── ground_truth.md              # Equation verification table
│   └── computational_verification.py
├── figures/
│   ├── exact_statmech.py            # ★ Core module: exact 2D stat-mech
│   ├── fig2_exact.py                # ★ Correct Figure 2 generator
│   ├── fig2_exact.pdf               # Exact Figure 2 (N=16)
│   ├── reproduce_figures.py         # Approximate figures (1D model)
│   ├── generalizations.py           # 6 generalization studies
│   ├── fig[2-6].pdf                 # Approximate paper figures
│   └── gen[1-6]_*.pdf               # Generalization figures
├── report/
│   ├── report.tex                   # Full LaTeX report
│   └── report.pdf                   # Compiled 16-page report
└── HANDOFF.md                       # This file
```
