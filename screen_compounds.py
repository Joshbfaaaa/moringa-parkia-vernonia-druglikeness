"""
RDKit-based drug-likeness / cheminformatics screening
Plant-derived antimicrobial compounds from Moringa oleifera, Parkia biglobosa,
and Vernonia amygdalina (companion analysis to Akpor et al. 2021, Sci. African).

USAGE (Google Colab or any machine with internet access):
    !pip install rdkit
    python screen_compounds.py compound_library.csv results_screening.csv

Requires: rdkit (pip install rdkit), pandas
"""

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def build_alert_catalog():
    """PAINS (Pan-Assay INterference compounds) + Brenk structural alert filters.
    Brenk's list flags reactive / toxicophoric / unstable substructures commonly
    used as a fast in-silico toxicity/liability screen ahead of assay work."""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog(params)


def lipinski_ro5(mw, logp, hbd, hba):
    violations = 0
    if mw > 500: violations += 1
    if logp > 5: violations += 1
    if hbd > 5: violations += 1
    if hba > 10: violations += 1
    return violations


def veber_rule(rotb, tpsa):
    # Veber et al. 2002: oral bioavailability heuristic
    return (rotb <= 10) and (tpsa <= 140)


def screen(input_csv, output_csv):
    lib = pd.read_csv(input_csv)
    catalog = build_alert_catalog()
    rows = []

    for _, r in lib.iterrows():
        mol = Chem.MolFromSmiles(r["smiles"])
        if mol is None:
            rows.append({**r.to_dict(), "ERROR": "invalid SMILES"})
            continue

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        rotb = Lipinski.NumRotatableBonds(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)
        arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

        ro5_violations = lipinski_ro5(mw, logp, hbd, hba)
        veber_pass = veber_rule(rotb, tpsa)

        matches = catalog.GetMatches(mol)
        alert_names = "; ".join(sorted({m.GetDescription() for m in matches})) or "none"

        rows.append({
            "name": r["name"],
            "plant": r["plant"],
            "class": r["phytochemical_class"],
            "MW": round(mw, 1),
            "LogP": round(logp, 2),
            "HBD": hbd,
            "HBA": hba,
            "TPSA": round(tpsa, 1),
            "RotatableBonds": rotb,
            "Rings": rings,
            "AromaticRings": arom_rings,
            "Lipinski_violations(of 4)": ro5_violations,
            "Passes_Ro5(<=1 violation)": ro5_violations <= 1,
            "Passes_Veber": veber_pass,
            "Drug-like(Ro5 & Veber)": (ro5_violations <= 1) and veber_pass,
            "Structural_alerts(PAINS/Brenk)": alert_names,
        })

    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)
    print(out.to_string(index=False))
    print(f"\nSaved: {output_csv}")
    print(f"\n{out['Drug-like(Ro5 & Veber)'].sum()} / {len(out)} compounds pass both Lipinski (<=1 violation) and Veber criteria.")
    return out


if __name__ == "__main__":
    in_csv = sys.argv[1] if len(sys.argv) > 1 else "compound_library.csv"
    out_csv = sys.argv[2] if len(sys.argv) > 2 else "results_screening.csv"
    screen(in_csv, out_csv)
