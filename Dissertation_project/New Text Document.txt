# ==============================
# Insect ligand–protein dataset
# Source: ChEMBL
# ==============================

from chembl_webresource_client.new_client import new_client
import pandas as pd

activity = new_client.activity
target = new_client.target
molecule = new_client.molecule

# Explicit insect species list
INSECT_SPECIES = [
    "Drosophila melanogaster",
    "Anopheles gambiae",
    "Aedes aegypti",
    "Bombyx mori",
    "Helicoverpa armigera",
    "Plutella xylostella",
    "Nilaparvata lugens"
]

records = []
missing_molecules = 0
total_activities = 0

print("Starting ChEMBL insect data extraction...\n")

for species in INSECT_SPECIES:
    print(f"Processing species: {species}")

    targets = target.filter(organism=species)
    print(f"  Targets found: {len(targets)}")

    for t in targets:
        target_id = t.get("target_chembl_id")
        target_name = t.get("pref_name")

        if not t.get("target_components"):
            continue

        for comp in t["target_components"]:
            accessions = comp.get("accession")
            if not accessions:
                continue

            uniprot_id = accessions[0]

            acts = activity.filter(
                target_chembl_id=target_id,
                standard_type__in=["IC50", "Ki", "EC50"],
                standard_relation="=",
                standard_units="nM"
            )

            for a in acts:
                total_activities += 1

                mol_id = a.get("molecule_chembl_id")
                value = a.get("standard_value")

                if mol_id is None or value is None:
                    continue

                mol = molecule.get(mol_id)
                if mol is None:
                    missing_molecules += 1
                    continue

                structures = mol.get("molecule_structures")
                if not structures:
                    missing_molecules += 1
                    continue

                smiles = structures.get("canonical_smiles")
                if smiles is None:
                    missing_molecules += 1
                    continue

                records.append({
                    "Ligand_SMILES": smiles,
                    "Target_UniProt": uniprot_id,
                    "Target_Name": target_name,
                    "Species": species,
                    "Activity_nM": float(value)
                })

print("\nExtraction finished.")
print("Total activity records scanned:", total_activities)
print("Missing / unusable molecules:", missing_molecules)
print("Collected ligand–target pairs:", len(records))

# Convert to DataFrame
df = pd.DataFrame(records)

if len(df) == 0:
    print("\n❌ No data collected.")
    print("➡️ Expand the species list or relax activity filters.")
else:
    # Remove duplicates
    df = df.drop_duplicates(subset=["Ligand_SMILES", "Target_UniProt"])

    # Optional: remove very weak binders (>50 µM)
    df = df[df["Activity_nM"] <= 50000]

    df.to_csv("insect_ligand_targets_raw.csv", index=False)
    print("\n✅ Saved: insect_ligand_targets_raw.csv")

    print("\nDataset summary:")
    print("Ligand–target pairs:", len(df))
    print("Unique ligands:", df["Ligand_SMILES"].nunique())
    print("Unique targets:", df["Target_UniProt"].nunique())
