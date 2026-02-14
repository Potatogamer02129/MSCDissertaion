import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.DataStructs import TanimotoSimilarity

# --- 1. PAGE CONFIGURATION ---
st.set_page_config(
    page_title="Target Prediction",
    page_icon="🎯",
    layout="wide"
)

# --- 2. THEME-AGNOSTIC CSS ---
st.markdown("""
<style>
    /* Adapts to Dark/Light mode automatically */
    .success-card {
        padding: 1rem;
        border-radius: 8px;
        background-color: #DCFCE7; 
        border: 1px solid #86EFAC;
        color: #14532D; 
        margin-bottom: 1rem;
    }
    div[data-testid="metric-container"] {
        border: 1px solid rgba(128, 128, 128, 0.2);
        padding: 10px;
        border-radius: 8px;
    }
    .stButton>button {
        width: 100%;
        border-radius: 8px;
        height: 3em;
        font-weight: bold;
    }
</style>
""", unsafe_allow_html=True)

# --- 3. DATA LOADING (Cached) ---
@st.cache_data(show_spinner=False)
def load_datasets():
    """Loads and prepares fingerprints for both datasets simultaneously."""
    # 1. ChEMBL Data
    try:
        chembl_df = pd.read_csv("insect_ligand_targets_raw.csv")
        chembl_df["Mol"] = chembl_df["Ligand_SMILES"].apply(Chem.MolFromSmiles)
        chembl_df = chembl_df[chembl_df["Mol"].notnull()]
        chembl_df["FP"] = chembl_df["Mol"].apply(
            lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
        )
    except FileNotFoundError:
        chembl_df = pd.DataFrame() # Fallback for testing UI

    # 2. Literature Data
    try:
        lit_df = pd.read_csv("literature_insect_bioactivity.csv")
        lit_df["Mol"] = lit_df["SMILES"].apply(Chem.MolFromSmiles)
        lit_df = lit_df[lit_df["Mol"].notnull()]
        lit_df["FP"] = lit_df["Mol"].apply(
            lambda m: AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
        )
    except FileNotFoundError:
        lit_df = pd.DataFrame()

    return chembl_df, lit_df

chembl_df, literature_df = load_datasets()

# --- 4. ALGORITHMIC LOGIC ---
def predict_targets_advanced(query_smiles, ref_df, cutoff=0.3):
    """Predicts targets using both Max and Mean Tanimoto similarity."""
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None:
        return None, "Invalid SMILES string."

    qfp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    results = []

    # Group by target to find cluster similarities
    for target, group in ref_df.groupby("Target_UniProt"):
        sims = [TanimotoSimilarity(qfp, rfp) for rfp in group["FP"]]
        valid_sims = [s for s in sims if s >= cutoff]

        if valid_sims:
            max_sim = max(valid_sims)
            mean_sim = sum(valid_sims) / len(valid_sims)
            
            # Upgraded Scoring: heavily weights the single closest neighbor (Max Sim) 
            # while still rewarding cluster density (Mean Sim & Count).
            score = (max_sim * 2) + (len(valid_sims) * mean_sim * 0.5)

            results.append({
                "Target Name": group["Target_Name"].iloc[0],
                "Species": group["Species"].iloc[0],
                "Active Ligands Matched": len(valid_sims),
                "Max Similarity": round(max_sim, 3),
                "Mean Similarity": round(mean_sim, 3),
                "Base Score": round(score, 3)
            })

    if not results:
        return pd.DataFrame(), None

    return pd.DataFrame(results).sort_values("Base Score", ascending=False).reset_index(drop=True), None

def dynamic_literature_bonus(query_smiles, lit_df):
    """Dynamically finds literature matches without hardcoding protein names."""
    if lit_df.empty: return 0.0, None
    
    mol = Chem.MolFromSmiles(query_smiles)
    if mol is None: return 0.0, None

    qfp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
    best_sim, best_row = 0.0, None

    for _, row in lit_df.iterrows():
        sim = TanimotoSimilarity(qfp, row["FP"])
        if sim > best_sim:
            best_sim = sim
            best_row = row

    # Define confidence tiers
    if best_sim >= 0.6:
        return 2.0, best_row  # High confidence boost
    elif best_sim >= 0.45:
        return 1.0, best_row  # Moderate confidence boost
    return 0.0, None

# --- 5. DASHBOARD UI ---

# Sidebar Controls
with st.sidebar:
    st.header("⚙️ Prediction Parameters")
    
    with st.expander("📌 Instructions", expanded=True):
        st.write("""
        1. Input a chemical SMILES string.
        2. Adjust the Tanimoto similarity threshold.
        3. Run prediction to mine ChEMBL and literature datasets.
        """)
        
    cutoff = st.slider(
        "Tanimoto Similarity Cutoff",
        min_value=0.2, max_value=0.8, value=0.35, step=0.05,
        help="Higher values require stricter chemical similarity to known ligands."
    )
    
    st.divider()
    st.caption("Module 1 | Chemoinformatics Engine")

# Main Content
st.title("🎯 Ligand-Based Target Prediction")
st.markdown("Predict probable insect protein targets for uncharacterized compounds using 2048-bit Morgan Fingerprints and structural similarity profiling.")

col1, col2 = st.columns([2, 1])

with col1:
    query_smiles = st.text_input(
        "Query Compound (SMILES)",
        placeholder="e.g. COC1=CC=CC=C1OC (Type or paste SMILES here)",
        help="Simplified Molecular-Input Line-Entry System"
    )
    run_btn = st.button("🚀 Execute Target Mining")

with col2:
    # 2D Molecule Visualization Box
    if query_smiles:
        mol = Chem.MolFromSmiles(query_smiles)
        if mol:
            # Generate the 2D image
            img = Draw.MolToImage(mol, size=(300, 150))
            st.image(img, caption="Query Structure Validated", use_container_width=True)
        else:
            st.error("Invalid SMILES structure.")

# --- 6. EXECUTION PHASE ---
if run_btn and query_smiles:
    if chembl_df.empty:
        st.error("Dataset Error: 'insect_ligand_targets_raw.csv' not found. Please ensure data is loaded.")
        st.stop()
        
    with st.status("Mining bioactivity databases...", expanded=True) as status:
        st.write("1️⃣ Encoding Query to 2048-bit Morgan Fingerprint...")
        
        # 1. Run ChEMBL Predictions
        st.write("2️⃣ Calculating Tanimoto similarities across ChEMBL dataset...")
        results_df, error = predict_targets_advanced(query_smiles, chembl_df, cutoff)
        
        if error:
            status.update(label="Execution Failed", state="error")
            st.error(error)
            st.stop()
            
        if results_df.empty:
            status.update(label="No Hits Found", state="error")
            st.warning(f"No targets found with a structural similarity above {cutoff}. Try lowering the cutoff.")
            st.stop()

        # 2. Run Dynamic Literature Bonus
        st.write("3️⃣ Cross-referencing against verified literature targets...")
        bonus_score, lit_hit = dynamic_literature_bonus(query_smiles, literature_df)
        
        # Apply the bonus dynamically to the Target Name matched in literature
        results_df["Literature Bonus"] = 0.0
        if lit_hit is not None:
            # Check if the literature target matches any of our ChEMBL predictions
            lit_target = lit_hit.get("Target_Name", "") # Safely get target name
            mask = results_df["Target Name"].str.contains(str(lit_target), case=False, na=False)
            results_df.loc[mask, "Literature Bonus"] = bonus_score

        # 3. Final Scoring & Formatting
        results_df["Final Predictive Score"] = results_df["Base Score"] + results_df["Literature Bonus"]
        
        # Confidence Labeling based on Max Similarity (much more accurate than Mean)
        def assign_confidence(row):
            if row["Max Similarity"] >= 0.6 or row["Literature Bonus"] > 0: return "🟢 High"
            elif row["Max Similarity"] >= 0.45: return "🟡 Medium"
            return "🔴 Low"
            
        results_df["Confidence"] = results_df.apply(assign_confidence, axis=1)
        
        # Clean up column order for display
        final_cols = ["Confidence", "Target Name", "Species", "Max Similarity", "Mean Similarity", "Active Ligands Matched", "Final Predictive Score"]
        display_df = results_df[final_cols].sort_values("Final Predictive Score", ascending=False)
        
        status.update(label="Target Mining Complete!", state="complete", expanded=False)

    # --- 7. RESULTS DASHBOARD ---
    st.divider()
    
    top_hit = display_df.iloc[0]
    
    # Top Metrics
    m1, m2, m3 = st.columns(3)
    m1.metric("Top Predicted Target", top_hit["Target Name"][:30])
    m2.metric("Max Similarity", f"{top_hit['Max Similarity']*100:.1f}%")
    m3.metric("Prediction Confidence", top_hit["Confidence"].replace("🟢 ", "").replace("🟡 ", "").replace("🔴 ", ""))
    
    # Data Tabs
    tab1, tab2 = st.tabs(["📊 Target Rankings", "📚 Literature Validation"])
    
    with tab1:
        st.markdown("#### Probable Protein Targets")
        st.dataframe(display_df, use_container_width=True, hide_index=True)
        
    with tab2:
        if lit_hit is not None:
            st.markdown(f"""
            <div class="success-card">
                <b>🔬 Verified Literature Match Found!</b><br>
                Your query compound exhibits structural similarity to a known active ligand.
            </div>
            """, unsafe_allow_html=True)
            
            # Recalculate exact sim for display
            q_mol = Chem.MolFromSmiles(query_smiles)
            q_fp = AllChem.GetMorganFingerprintAsBitVect(q_mol, 2, 2048)
            exact_sim = TanimotoSimilarity(q_fp, lit_hit['FP'])
            
            st.write(f"**Matched Compound:** {lit_hit.get('Ligand_Name', 'Unknown')}")
            st.write(f"**Target Pathway:** {lit_hit.get('Target_Name', 'Unknown')}")
            st.write(f"**Structural Similarity:** {exact_sim:.3f}")
            st.write(f"**Evidence:** {lit_hit.get('Evidence_Type', 'Literature')}")
            st.caption(f"Citation: {lit_hit.get('Reference', 'N/A')}")
        else:
            st.info("No high-confidence matches found in the literature-curated dataset. Predictions rely entirely on ChEMBL data.")
