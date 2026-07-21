"""
Generate summary figures from results_screening.csv
(Companion to screen_compounds.py — Moringa/Parkia/Vernonia drug-likeness screen)

USAGE:
    python3 make_figures.py results_screening.csv
Outputs (saved in same folder):
    fig1_druglikeness_summary.png
    fig2_lipinski_space.png
    fig3_logp_by_compound.png
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt

def main(csv_path="results_screening.csv"):
    df = pd.read_csv(csv_path)
    df = df.dropna(subset=["MW"])  # drop any invalid-SMILES rows

    # ---- Fig 1: Pass/fail summary bar ----
    fig, ax = plt.subplots(figsize=(5, 4))
    counts = df["Drug-like(Ro5 & Veber)"].value_counts()
    labels = ["Pass" if v else "Fail" for v in counts.index]
    colors = ["#4C956C" if l == "Pass" else "#D64550" for l in labels]
    ax.bar(labels, counts.values, color=colors)
    for i, v in enumerate(counts.values):
        ax.text(i, v + 0.1, str(v), ha="center", fontweight="bold")
    ax.set_ylabel("Number of compounds")
    ax.set_title("Drug-likeness screen: Lipinski (\u22641 violation) + Veber")
    plt.tight_layout()
    plt.savefig("fig1_druglikeness_summary.png", dpi=200)
    plt.close()

    # ---- Fig 2: Chemical space scatter (MW vs TPSA), colored by pass/fail ----
    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    for passed, group in df.groupby("Drug-like(Ro5 & Veber)"):
        color = "#4C956C" if passed else "#D64550"
        label = "Drug-like" if passed else "Fails Ro5/Veber"
        ax.scatter(group["MW"], group["TPSA"], s=90, color=color, edgecolor="black",
                   linewidth=0.5, label=label, zorder=3)

    # Detect points that sit on top of one another (rounded MW/TPSA match) and
    # fan their labels out in a wide arc with thin leader lines so text never overlaps.
    # Also treat the dense mid-cluster (Quercetin/Kaempferol/Luteolin/Rutin/Catechin/
    # Epicatechin, all MW 280-300, TPSA 110-135) as one group so their labels spread
    # further apart rather than only de-duplicating exact ties.
    df["_cluster"] = "solo"
    mid_mask = df["MW"].between(280, 320) & df["TPSA"].between(105, 140)
    df.loc[mid_mask, "_cluster"] = "mid"
    sterol_mask = df["MW"].between(400, 430) & df["TPSA"].between(15, 25)
    df.loc[sterol_mask, "_cluster"] = "sterol"

    wide_offsets = [(30, 55), (85, 15), (95, -35), (30, -70), (-90, -55), (-95, 20)]
    sterol_offsets = [(-90, 30), (20, 55), (60, -10)]

    for _, r in df[df["_cluster"] == "solo"].iterrows():
        ax.annotate(
            r["name"], (r["MW"], r["TPSA"]), xytext=(20, 14),
            textcoords="offset points", fontsize=7.5, ha="left",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.6, shrinkA=3, shrinkB=3),
        )

    mid_group = df[df["_cluster"] == "mid"].sort_values("TPSA")
    for i, (_, r) in enumerate(mid_group.iterrows()):
        dx, dy = wide_offsets[i % len(wide_offsets)]
        ax.annotate(
            r["name"], (r["MW"], r["TPSA"]),
            xytext=(dx, dy), textcoords="offset points",
            fontsize=7.5, ha="left" if dx >= 0 else "right",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.6, shrinkA=3, shrinkB=3),
        )

    sterol_group = df[df["_cluster"] == "sterol"].sort_values("MW")
    for i, (_, r) in enumerate(sterol_group.iterrows()):
        dx, dy = sterol_offsets[i % len(sterol_offsets)]
        ax.annotate(
            r["name"], (r["MW"], r["TPSA"]),
            xytext=(dx, dy), textcoords="offset points",
            fontsize=7.5, ha="left" if dx >= 0 else "right",
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.6, shrinkA=3, shrinkB=3),
        )

    ax.axvline(500, color="gray", linestyle="--", linewidth=1, label="MW = 500 (Ro5 limit)")
    ax.axhline(140, color="gray", linestyle=":", linewidth=1, label="TPSA = 140 (Veber limit)")
    ax.set_xlabel("Molecular Weight (Da)")
    ax.set_ylabel("TPSA (\u00c5\u00b2)")
    ax.set_title("Chemical space: MW vs TPSA")
    ax.set_xlim(df["MW"].min() - 110, df["MW"].max() + 100)
    ax.set_ylim(df["TPSA"].min() - 40, df["TPSA"].max() + 40)
    ax.legend(fontsize=8, loc="upper left")
    plt.tight_layout()
    plt.savefig("fig2_lipinski_space.png", dpi=200)
    plt.close()

    # ---- Fig 3: LogP per compound, sorted ----
    fig, ax = plt.subplots(figsize=(8, 5.5))
    d = df.sort_values("LogP")
    colors = ["#D64550" if lp > 5 else "#4C956C" for lp in d["LogP"]]
    ax.barh(d["name"], d["LogP"], color=colors)
    ax.axvline(5, color="black", linestyle="--", linewidth=1, label="LogP = 5 (Ro5 limit)")
    ax.set_xlabel("Calculated LogP (Crippen)")
    ax.set_title("Lipophilicity across the screened phytochemical library")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig("fig3_logp_by_compound.png", dpi=200)
    plt.close()

    print("Saved: fig1_druglikeness_summary.png, fig2_lipinski_space.png, fig3_logp_by_compound.png")

if __name__ == "__main__":
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "results_screening.csv"
    main(csv_path)
