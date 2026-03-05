#!/usr/bin/env python3
"""
Plot FeO potential energy curves from compute_feo_pec.py output.

Usage:
  python plot_feo_pec.py                          # NEVPT2 energies
  python plot_feo_pec.py --level casscf           # CASSCF energies
  python plot_feo_pec.py --max-energy 35000       # truncate at 35000 cm^-1
"""

import numpy as np
import json
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Hartree to cm^-1
HA_TO_CM = 219474.63

# C2v -> Lambda-S label mapping (approximate — based on root ordering)
# The actual assignment requires CI vector analysis; this gives conventional labels.
LAMBDA_LABELS = {
    'A1': ['Sigma+', 'Delta', 'Gamma'],
    'A2': ['Sigma-', 'Delta', 'Gamma'],
    'B1': ['Pi', 'Phi'],
    'B2': ['Pi', 'Phi'],
}

MULT_COLORS = {
    3: plt.cm.Blues,
    5: plt.cm.Reds,
    7: plt.cm.Greens,
}

MULT_NAMES = {3: 'Triplet', 5: 'Quintet', 7: 'Septet'}


def load_data(result_file):
    with open(result_file) as f:
        data = json.load(f)
    return data['config'], data['results']


def extract_curves(results, level='nevpt2'):
    """
    Extract PEC data as a list of dicts:
      {mult, irrep, root, r_values, energies, label}
    """
    # Collect all bond lengths
    r_keys = sorted(results.keys(), key=float)

    # Find all state blocks
    sample = results[r_keys[0]]
    state_keys = sorted(sample.keys())

    curves = []
    for sk in state_keys:
        mult, irrep = sk.split('_')
        mult = int(mult)
        nroots = len(sample[sk][level])

        for iroot in range(nroots):
            rs = []
            es = []
            for rk in r_keys:
                block = results[rk].get(sk)
                if block is None:
                    continue
                e = block[level][iroot]
                if e is not None:
                    rs.append(float(rk))
                    es.append(e)

            if len(rs) < 3:
                continue

            # Lambda label
            lambda_labels = LAMBDA_LABELS.get(irrep, [])
            if iroot < len(lambda_labels):
                lbl = lambda_labels[iroot]
            else:
                lbl = f'root{iroot}'

            sup = {3: '3', 5: '5', 7: '7'}.get(mult, str(mult))
            label = f'$^{sup}${lbl}'

            curves.append({
                'mult': mult,
                'irrep': irrep,
                'root': iroot,
                'r': np.array(rs),
                'e': np.array(es),
                'label': label,
                'state_key': sk,
            })

    return curves


def identify_degenerate_pairs(curves, tol_cm=50):
    """
    Identify degenerate Pi, Delta, Phi pairs (B1/B2 or A1/A2)
    and mark duplicates for removal from legend.
    """
    pairs = {}
    for c in curves:
        key = (c['mult'], c['root'], c['label'])
        pairs.setdefault(key, []).append(c)

    dedup = set()
    for key, group in pairs.items():
        if len(group) == 2:
            # Check if they're close in energy at equilibrium
            r_common = np.intersect1d(
                np.round(group[0]['r'], 4),
                np.round(group[1]['r'], 4)
            )
            if len(r_common) > 0:
                r0 = r_common[len(r_common)//2]
                idx0 = np.argmin(np.abs(group[0]['r'] - r0))
                idx1 = np.argmin(np.abs(group[1]['r'] - r0))
                de = abs(group[0]['e'][idx0] - group[1]['e'][idx1]) * HA_TO_CM
                if de < tol_cm:
                    # Mark second one as duplicate
                    dedup.add(id(group[1]))
    return dedup


def find_spectroscopic_constants(r_ang, e_hartree):
    """
    Fit Morse potential near minimum to extract re, De, omega_e.
    Returns (re_ang, De_cm, omega_e_cm) or None if no minimum.
    """
    if len(r_ang) < 5:
        return None

    imin = np.argmin(e_hartree)
    if imin == 0 or imin == len(e_hartree) - 1:
        return None  # no interior minimum

    # Fit 4th-order polynomial around minimum
    lo = max(0, imin - 4)
    hi = min(len(r_ang), imin + 5)
    r_fit = r_ang[lo:hi]
    e_fit = e_hartree[lo:hi]

    if len(r_fit) < 5:
        return None

    coeffs = np.polyfit(r_fit, e_fit, 4)
    p = np.poly1d(coeffs)
    dp = p.deriv()

    # Find minimum from derivative roots
    roots = dp.roots
    real_roots = roots[np.isreal(roots)].real
    real_roots = real_roots[(real_roots > r_fit[0]) & (real_roots < r_fit[-1])]
    if len(real_roots) == 0:
        return None

    re = real_roots[np.argmin(np.abs(real_roots - r_ang[imin]))]
    e_min = p(re)

    # omega_e from second derivative: omega_e = sqrt(k/mu) where k = d2E/dr2
    d2p = dp.deriv()
    k = d2p(re)  # in Hartree / Ang^2

    # Reduced mass of FeO in atomic units
    m_Fe = 55.845  # amu
    m_O = 15.999
    mu_amu = m_Fe * m_O / (m_Fe + m_O)
    mu_au = mu_amu * 1822.888  # amu -> electron masses

    # Convert k to atomic units: Hartree / bohr^2
    k_au = k * (0.529177**2)  # Ang^2 -> bohr^2

    if k_au <= 0:
        return None

    omega_au = np.sqrt(k_au / mu_au)
    omega_cm = omega_au * HA_TO_CM

    # Dissociation energy (crude: E(r_max) - E(re))
    De = (e_hartree[-1] - e_min) * HA_TO_CM

    return re, De, omega_cm


def plot_curves(curves, config, level, max_energy_cm=None, outfile=None):
    """Create publication-quality PEC plot."""
    if outfile is None:
        outfile = Path('feo_pec_results') / f'feo_pec_{level}.png'

    # Find global minimum for energy reference
    e_min_global = min(np.min(c['e']) for c in curves)

    dedup = identify_degenerate_pairs(curves)

    fig, ax = plt.subplots(figsize=(12, 8))

    legend_handles = []
    legend_labels = []
    plotted_labels = set()

    for c in curves:
        mult = c['mult']
        cmap = MULT_COLORS.get(mult, plt.cm.Greys)

        # Color: spread roots across colormap
        nroots_this = sum(1 for cc in curves
                          if cc['mult'] == mult and cc['irrep'] == c['irrep'])
        frac = 0.35 + 0.55 * c['root'] / max(nroots_this - 1, 1)
        color = cmap(frac)

        # Linestyle by irrep type
        if c['irrep'] in ('A1', 'A2'):
            ls = '-'
        else:
            ls = '--'

        e_rel = (c['e'] - e_min_global) * HA_TO_CM

        if max_energy_cm and np.min(e_rel) > max_energy_cm:
            continue

        ax.plot(c['r'], e_rel, color=color, ls=ls, lw=1.5, alpha=0.85)

        # Label at equilibrium (or leftmost point)
        imin = np.argmin(c['e'])
        if id(c) not in dedup and c['label'] not in plotted_labels:
            ax.annotate(
                c['label'], (c['r'][imin], e_rel[imin]),
                fontsize=7, ha='left', va='bottom',
                xytext=(4, 4), textcoords='offset points',
                alpha=0.8,
            )
            plotted_labels.add(c['label'])

    # Custom legend for multiplicities
    for mult, name in sorted(MULT_NAMES.items()):
        cmap = MULT_COLORS[mult]
        legend_handles.append(Line2D([0], [0], color=cmap(0.6), lw=2))
        legend_labels.append(f'{name} (2S+1={mult})')

    legend_handles.append(Line2D([0], [0], color='gray', ls='-', lw=1.5))
    legend_labels.append(r'$\Sigma, \Delta, \Gamma$ (A1/A2)')
    legend_handles.append(Line2D([0], [0], color='gray', ls='--', lw=1.5))
    legend_labels.append(r'$\Pi, \Phi$ (B1/B2)')

    ax.legend(legend_handles, legend_labels, loc='upper right', fontsize=9)

    if max_energy_cm:
        ax.set_ylim(-500, max_energy_cm)

    ax.set_xlabel(r'$r_{\mathrm{Fe-O}}$ ($\mathrm{\AA}$)', fontsize=13)
    ax.set_ylabel(r'Energy (cm$^{-1}$)', fontsize=13)
    ax.set_title(
        f'FeO Potential Energy Curves — '
        f'SA-CASSCF({config["nelecas"]},{config["ncas"]})+{level.upper()}'
        f'/{config["basis"]} (X2C)',
        fontsize=12,
    )
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(outfile, dpi=200)
    print(f"Plot saved: {outfile}")

    # Also save PDF
    pdf_file = outfile.with_suffix('.pdf')
    fig.savefig(pdf_file)
    print(f"Plot saved: {pdf_file}")
    plt.close(fig)


def print_spectroscopic_table(curves, e_min_global):
    """Print table of spectroscopic constants for bound states."""
    print("\n" + "="*75)
    print("  Spectroscopic Constants (NEVPT2)")
    print("="*75)
    print(f"  {'State':<16} {'re (A)':>8} {'Te (cm-1)':>12} "
          f"{'omega_e (cm-1)':>15} {'De (cm-1)':>12}")
    print("-"*75)

    dedup = identify_degenerate_pairs(curves)
    state_data = []

    for c in curves:
        if id(c) in dedup:
            continue
        sc = find_spectroscopic_constants(c['r'], c['e'])
        if sc is None:
            continue
        re, De, omega_e = sc
        Te = (np.min(c['e']) - e_min_global) * HA_TO_CM
        state_data.append((Te, c['label'], re, omega_e, De))

    state_data.sort()
    for Te, label, re, omega_e, De in state_data[:25]:
        print(f"  {label:<16} {re:8.4f} {Te:12.0f} {omega_e:15.0f} {De:12.0f}")

    print("="*75)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--level', default='nevpt2', choices=['casscf', 'nevpt2'])
    parser.add_argument('--max-energy', type=float, default=45000,
                        help='Max energy in cm^-1 for plot (default: 45000)')
    parser.add_argument('--data', default='feo_pec_results/pec_data.json')
    args = parser.parse_args()

    config, results = load_data(args.data)
    curves = extract_curves(results, level=args.level)

    if not curves:
        print("No data to plot!")
        raise SystemExit(1)

    e_min_global = min(np.min(c['e']) for c in curves)

    print(f"Loaded {len(curves)} curves, {len(results)} bond lengths")
    print(f"Level: {args.level.upper()}")
    print(f"Basis: {config['basis']}")

    plot_curves(curves, config, args.level, max_energy_cm=args.max_energy)
    print_spectroscopic_table(curves, e_min_global)


if __name__ == '__main__':
    main()
