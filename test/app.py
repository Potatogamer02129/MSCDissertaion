import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity

# =========================
# App config
# =========================
st.set_page_config(
    page_title="Insect Target Prediction",
    layout="centered"
)

st.title("🦟 Insect Protein Target Prediction")
st.write(
    "Ligand-based prediction of probable insect protein targets "
    "using experimentally validated ChEMBL insect bioactivity data."
)

# =========================
# Load + clean + fingerprint (ONCE)
# =========================
@st.cache_resource
def load_and_prepare_data():
    df = pd.read_csv("insect_ligand_targets_raw.csv")

    # --- Keep only required columns ---
    keep_cols = [
        "Ligand_SMILES",
        "Target_CHEMBL",
        "Target_Name",
        "Species",
        "Activity_Type",
        "Activity_Value"
    ]
    df = df[keep_cols]

    # --- Validate SMILES ---
    def mol_from_smiles(s):
        try:
            return Chem.MolFromSmiles(s)
        except:
            return None

    df["Mol"] = df["Ligand_SMILES"].apply(mol_from_smiles)
    df = df[df["Mol"].notnull()]

    # --- Deduplicate ---
    df = df.drop_duplicates(
        subset=["Ligand_SMILES", "Target_CHEMBL"]
    )

    # --- Fingerprints ---
    df["FP"] = df["Mol"].apply(
        lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
    )

    return df.reset_index(drop=True)

df = load_and_prepare_data()

st.success(
    f"Dataset ready: {df['Ligand_SMILES'].nunique()} ligands | "
    f"{df['Target_CHEMBL'].nunique()} insect targets"
)

# =========================
# Prediction function
# =========================
def predict_targets(query_smiles, cutoff):
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        return None, "Invalid SMILES"

    qfp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    results = []

    for tid, group in df.groupby("Target_CHEMBL"):
        sims = [
            TanimotoSimilarity(qfp, fp)
            for fp in group["FP"]
            if TanimotoSimilarity(qfp, fp) >= cutoff
        ]

        if sims:
            mean_sim = sum(sims) / len(sims)
            score = mean_sim * len(sims)

            results.append({
                "Target name": group["Target_Name"].iloc[0],
                "Species": group["Species"].iloc[0],
                "Matched ligands": len(sims),
                "Mean similarity": round(mean_sim, 3),
                "Score": round(score, 3)
            })

    if not results:
        return pd.DataFrame(), None

    return (
        pd.DataFrame(results)
        .sort_values("Score", ascending=False)
        .reset_index(drop=True),
        None
    )

# =========================
# UI
# =========================
st.subheader("🔬 Query compound")

query_smiles = st.text_input(
    "Enter SMILES",
    placeholder="e.g. COC1=CC=CC=C1OC"
)

cutoff = st.slider(
    "Tanimoto similarity cutoff",
    min_value=0.2,
    max_value=0.6,
    value=0.3,
    step=0.05
)

if st.button("Predict insect targets"):
    if not query_smiles.strip():
        st.warning("Please enter a SMILES string.")
    else:
        results, error = predict_targets(query_smiles, cutoff)

        if error:
            st.error(error)
        elif results.empty:
            st.info("No insect protein targets found above the cutoff.")
        else:
            st.subheader("📊 Predicted insect targets")
            st.dataframe(results)

            st.caption(
                "⚠️ Predictions are based on chemical similarity to known insect ligands. "
                "They do not imply binding affinity or biological activity."
            )
