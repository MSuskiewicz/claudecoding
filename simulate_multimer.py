#!/usr/bin/env python3
"""
simulate_multimer.py

Compute and plot the equilibrium distribution of oligomeric species
(1-mer through N-mer) for a filamentous multimer with a constant
step-wise dissociation constant (Kd) between each protomer addition.

Model
-----
  A1 + A1  <-->  A2    Kd = [A1]^2 / [A2]
  A2 + A1  <-->  A3    Kd = [A2][A1] / [A3]
  ...
  A(n-1) + A1  <-->  An

This gives  [An] = [A1]^n / Kd^(n-1)

Conservation: C_total = sum_{n=1}^{N} n * [An]

[A1] is found by numerical root-finding; all other species follow analytically.

Usage
-----
  python simulate_multimer.py --Kd 1.0 --Ctotal 10.0 --Nmax 12

Arguments
---------
  --Kd      Dissociation constant in µM  (default: 1.0)
  --Ctotal  Total protein concentration in µM  (default: 10.0)
  --Nmax    Largest oligomer to include  (default: 12)
  --out     Output image filename  (default: multimer_distribution.png)
  --yaxis   Y-axis quantity: 'fraction' (fraction of total protein mass)
            or 'conc' (molar concentration in µM)  (default: fraction)
"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


def species_concentrations(m1, Kd, Nmax):
    """
    Return list of concentrations [A1], [A2], ..., [ANmax] in µM.

    [An] = m1^n / Kd^(n-1)
    """
    return [m1 ** n / Kd ** (n - 1) for n in range(1, Nmax + 1)]


def total_monomer_equiv(m1, Kd, Nmax):
    """
    Sum_{n=1}^{Nmax} n * [An]  --  monomer-equivalent total concentration.
    """
    cn = species_concentrations(m1, Kd, Nmax)
    return sum((n + 1) * c for n, c in enumerate(cn))


def solve_monomer(Ctotal, Kd, Nmax):
    """
    Find the free monomer concentration [A1] that satisfies the conservation
    equation by bisection.
    """
    # Lower bound: m1 -> 0  =>  total -> 0
    # Upper bound: all protein is monomeric  =>  m1 = Ctotal
    lo, hi = 1e-12, Ctotal

    # Sanity: at hi the total_monomer_equiv >= Ctotal (monotone function)
    if total_monomer_equiv(hi, Kd, Nmax) < Ctotal:
        sys.exit(
            "Conservation equation cannot be satisfied.  "
            "Try a smaller Nmax or a larger Kd."
        )

    m1 = brentq(lambda m: total_monomer_equiv(m, Kd, Nmax) - Ctotal,
                lo, hi, xtol=1e-14, rtol=1e-12)
    return m1


def plot_distribution(concentrations, Ctotal, Kd, Nmax, yaxis, outfile):
    n_values = list(range(1, Nmax + 1))

    if yaxis == "fraction":
        # Fraction of total protein mass in each oligomeric state
        monomer_equiv = [n * c for n, c in zip(n_values, concentrations)]
        total = sum(monomer_equiv)
        y = [v / total for v in monomer_equiv]
        ylabel = "Fraction of total protein"
        title_extra = ""
    else:
        y = concentrations
        ylabel = "Concentration (µM)"
        title_extra = ""

    fig, ax = plt.subplots(figsize=(max(6, Nmax * 0.7), 5))
    bars = ax.bar(n_values, y, color="steelblue", edgecolor="white", linewidth=0.5)

    ax.set_xlabel("Oligomeric state (n-mer)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(
        f"Multimer distribution\n"
        f"Kd = {Kd:.3g} µM   |   C_total = {Ctotal:.3g} µM   |   max n-mer = {Nmax}",
        fontsize=11,
    )
    ax.set_xticks(n_values)
    ax.set_xticklabels([f"{n}" for n in n_values])

    # Annotate bars with values > 1 % (or > 1 nM for conc mode)
    threshold = 0.01 if yaxis == "fraction" else 0.001
    for bar, val in zip(bars, y):
        if val >= threshold:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(y) * 0.01,
                f"{val:.2f}" if yaxis == "fraction" else f"{val:.3g}",
                ha="center", va="bottom", fontsize=8,
            )

    ax.set_ylim(0, max(y) * 1.15)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150)
    print(f"Plot saved to {outfile}")
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Equilibrium distribution of filamentous multimers.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--Kd",     type=float, default=1.0,
                        help="Step-wise dissociation constant in µM (default: 1.0)")
    parser.add_argument("--Ctotal", type=float, default=10.0,
                        help="Total protein concentration in µM (default: 10.0)")
    parser.add_argument("--Nmax",   type=int,   default=12,
                        help="Largest oligomer to include (default: 12)")
    parser.add_argument("--out",    default="multimer_distribution.png",
                        help="Output image file (default: multimer_distribution.png)")
    parser.add_argument("--yaxis",  choices=["fraction", "conc"], default="fraction",
                        help="Y-axis: 'fraction' of total protein or 'conc' in µM (default: fraction)")
    args = parser.parse_args()

    if args.Kd <= 0:
        sys.exit("--Kd must be positive.")
    if args.Ctotal <= 0:
        sys.exit("--Ctotal must be positive.")
    if args.Nmax < 2:
        sys.exit("--Nmax must be >= 2.")

    print(f"Parameters: Kd = {args.Kd} µM, C_total = {args.Ctotal} µM, Nmax = {args.Nmax}")

    m1 = solve_monomer(args.Ctotal, args.Kd, args.Nmax)
    concentrations = species_concentrations(m1, args.Kd, args.Nmax)

    print(f"\n{'n-mer':>6}  {'[An] (µM)':>14}  {'mass fraction':>14}")
    print("-" * 40)
    total_mass = sum(n * c for n, c in enumerate(concentrations, 1))
    for n, c in enumerate(concentrations, 1):
        frac = n * c / total_mass
        print(f"{n:>6}  {c:>14.4g}  {frac:>14.4f}")
    print("-" * 40)
    print(f"{'Total':>6}  {sum(concentrations):>14.4g}  {'1.0000':>14}")

    plot_distribution(concentrations, args.Ctotal, args.Kd, args.Nmax,
                      args.yaxis, args.out)


if __name__ == "__main__":
    main()
