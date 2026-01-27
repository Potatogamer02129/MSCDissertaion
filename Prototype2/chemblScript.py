from chembl_webresource_client.new_client import new_client
import pandas as pd
import os
import time

print("Starting ChEMBL Insecta extraction (robust mode)...")

# Ensure data directory exists
os.makedirs("data", exist_ok=True)

activity = new_client.activity
target = new_client.target
molecule = new_client.molecule

records = []

# Correct taxonomic filter
insect_targets = target.filter(target_tax_id=50557)
total_targets = len(insect_targets)

print(f"Total insect targets found: {total_targets}")

for idx, t in enumerate(insect_targets, start=1):
    tid = t.get("target_chembl_id")
    if not tid:
        continue

    tname = t.get("pref_name", "")
    organism = t.get("organism", "")

    try:
        acts = activity.filter(
            target_chembl_id=tid,
            standard_type__in=["IC50", "EC50", "Ki", "Kd", "Inhibition"]
        )
    except Exception:
        continue

    for a in acts:
        mol_id = a.get("molecule_chembl_id")
        if not mol_id:
            continue

        try:
            mol = molecule.get(mol_id)
            if mol is None:
                continue

            structs = mol.get("molecule_structures")
            if not structs:
                continue

            smiles = structs.get("canonical_smiles")
            if not smiles:
                continue

        except Exception:
            continue

        records.append({
            "Ligand_SMILES": smiles,
            "Target_CHEMBL": tid,
            "Target_Name": tname,
            "Species": organism,
            "Activity_Type": a.get("standard_type"),
            "Activity_Value": a.get("standard_value"),
            "Activity_Unit": a.get("standard_units")
        })

        # Incremental save every 2000 records
        if len(records) % 2000 == 0:
            pd.DataFrame(records).drop_duplicates().to_csv(
                "data/insect_chembl_raw_partial.csv",
                index=False
            )

    # Progress indicator every 500 targets
    if idx % 500 == 0:
        print(f"Processed {idx}/{total_targets} targets | records: {len(records)}")

# Final save
df = pd.DataFrame(records).drop_duplicates()
df.to_csv("data/insect_chembl_raw.csv", index=False)

print("Extraction finished.")
print(f"Final ligand–target pairs: {len(df)}")
print("Saved: data/insect_chembl_raw.csv")
