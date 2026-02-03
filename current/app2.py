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
# Load ChEMBL insect dataset
# =========================
@st.cache_data
def load_chembl_dataset():
    df = pd.read_csv("insect_ligand_targets_raw.csv")
    df["Mol"] = df["Ligand_SMILES"].apply(Chem.MolFromSmiles)
    df = df[df["Mol"].notnull()]
    df["FP"] = df["Mol"].apply(
        lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
    )
    return df

# =========================
# Load literature AChE dataset
# =========================
@st.cache_data
def load_literature_dataset():
    df = pd.read_csv("literature_compounds.csv")
    df["Mol"] = df["SMILES"].apply(Chem.MolFromSmiles)
    df = df[df["Mol"].notnull()]
    df["FP"] = df["Mol"].apply(
        lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
    )
    return df

chembl_df = load_chembl_dataset()
literature_df = load_literature_dataset()

st.success(
    f"Loaded {chembl_df['Ligand_SMILES'].nunique()} ChEMBL ligands "
    f"and {literature_df.shape[0]} literature-curated compounds"
)

# =========================
# Core target prediction
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
                "ChEMBL score": round(score, 3)
            })

    if not results:
        return pd.DataFrame(), None

    return (
        pd.DataFrame(results)
        .sort_values("ChEMBL score", ascending=False)
        .reset_index(drop=True),
        None
    )

# =========================
# Literature AChE bonus
# =========================
def ache_literature_bonus(query_smiles, literature_df):
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        return 0.0, None

    qfp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    best_sim = 0
    best_row = None

    for _, row in literature_df.iterrows():
        sim = TanimotoSimilarity(qfp, row["FP"])
        if sim > best_sim:
            best_sim = sim
            best_row = row

    if best_sim >= 0.6:
        return 0.5, best_row
    elif best_sim >= 0.45:
        return 0.3, best_row
    else:
        return 0.0, None

# =========================
# Confidence labeling
# =========================
def confidence_label(row):
    if row["Matched ligands"] >= 5 and row["Mean similarity"] >= 0.35:
        return "Medium"
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
        results, error = predict_targets(query_smiles, chembl_df, cutoff)

        if error:
            st.error(error)

        elif results.empty:
            st.info(
                "No insect protein targets found above the selected similarity cutoff."
            )

        else:
            # Confidence
            results["Confidence"] = results.apply(confidence_label, axis=1)

            # Literature bonus
            bonus, lit_hit = ache_literature_bonus(query_smiles, literature_df)
            results["Literature bonus"] = 0.0

            ache_mask = results["Target name"].str.contains(
                "acetylcholinesterase", case=False
            )
            results.loc[ache_mask, "Literature bonus"] = bonus

            results["Final score"] = (
                results["ChEMBL score"] + results["Literature bonus"]
            )

            st.subheader("📊 Predicted insect protein targets")
            st.dataframe(
                results.sort_values("Final score", ascending=False)
            )

            # Show literature evidence separately
            if lit_hit is not None:
                st.subheader("📚 Supporting literature evidence (AChE)")
                st.markdown(
                    f"""
                    **Matched compound:** {lit_hit['Ligand_Name']}  
                    **Similarity:** {round(TanimotoSimilarity(
                        AllChem.GetMorganFingerprintAsBitVect(
                            Chem.MolFromSmiles(query_smiles), 2, 2048
                        ),
                        lit_hit['FP']
                    ), 3)}  
                    **Evidence:** {lit_hit['Evidence_Type']}  
                    **Citation:** {lit_hit['Cited_Source']}
                    """
                )

            st.caption(
                "⚠️ Predictions are based on chemical similarity and literature support. "
                "They do not imply confirmed binding or bioactivity."
            )

