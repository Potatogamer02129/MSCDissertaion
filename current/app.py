import streamlit as st
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity

# =========================
# App configuration
# =========================
st.set_page_config(
    page_title="Insect Target Prediction",
    layout="centered"
)

st.title("🦟 Insect Protein Target Prediction")
st.write(
    "Ligand-based prediction of **probable insect protein targets** "
    "using chemical similarity to known insect bioactive compounds."
)

# =========================
# Load dataset ONCE
# =========================
@st.cache_data
def load_dataset():
    df = pd.read_csv("insect_ligand_targets_raw.csv")
    df["Mol"] = df["Ligand_SMILES"].apply(Chem.MolFromSmiles)
    df = df[df["Mol"].notnull()]
    df["FP"] = df["Mol"].apply(
        lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
    )
    return df

df = load_dataset()

st.success(
    f"Dataset loaded: {df['Ligand_SMILES'].nunique()} ligands, "
    f"{df['Target_UniProt'].nunique()} insect protein targets"
)

# =========================
# Core prediction function
# =========================
def predict_targets(query_smiles, reference_df, cutoff=0.3):
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        return None, "Invalid SMILES string."

    qfp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    results = []

    for target, group in reference_df.groupby("Target_UniProt"):
        sims = []

        for rfp in group["FP"]:
            sim = TanimotoSimilarity(qfp, rfp)
            if sim >= cutoff:
                sims.append(sim)

        if sims:
            mean_sim = sum(sims) / len(sims)
            score = len(sims) * mean_sim

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
# Confidence labeling
# =========================
def confidence_label(row):
    if row["Matched ligands"] >= 5 and row["Mean similarity"] >= 0.35:
        return "Medium"
    else:
        return "Low"

# =========================
# User input
# =========================
st.subheader("🔬 Query compound")

query_smiles = st.text_input(
    "Enter compound SMILES",
    placeholder="e.g. COC1=CC=CC=C1OC"
)

cutoff = st.slider(
    "Tanimoto similarity cutoff",
    min_value=0.2,
    max_value=0.6,
    value=0.3,
    step=0.05
)

# =========================
# Run prediction
# =========================
if st.button("Predict insect targets"):
    if not query_smiles.strip():
        st.warning("Please enter a SMILES string.")
    else:
        results, error = predict_targets(query_smiles, df, cutoff)

        if error:
            st.error(error)
        elif results.empty:
            st.info(
                "No insect protein targets found above the selected similarity cutoff."
            )
        else:
            results["Confidence"] = results.apply(confidence_label, axis=1)

            st.subheader("📊 Predicted insect protein targets")
            st.dataframe(results)

            st.caption(
                "⚠️ Predictions are based on chemical similarity to known insect ligands. "
                "They do not imply binding affinity or biological activity."
            )

