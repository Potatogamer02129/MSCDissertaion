import streamlit as st

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="Dissertation Portal",
    page_icon="🎓",
    layout="wide"
)

# --- CUSTOM CSS FOR "HERO" SECTION ---
st.markdown("""
<style>
    /* Hero Title */
    .hero-title {
        font-size: 3rem;
        font-weight: 800;
        color: #0F172A;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    /* Subtitle */
    .hero-subtitle {
        font-size: 1.2rem;
        color: #64748B;
        text-align: center;
        margin-bottom: 2rem;
    }
    /* Card Container */
    .card {
        background-color: #FFFFFF;
        padding: 2rem;
        border-radius: 12px;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);
        border: 1px solid #E2E8F0;
        text-align: center;
        height: 100%;
        transition: transform 0.2s;
    }
    .card:hover {
        transform: translateY(-5px);
        box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
    }
    /* Card Icon */
    .icon {
        font-size: 3rem;
        margin-bottom: 1rem;
    }
    /* Button styling is handled by Streamlit, but we wrap it for alignment */
</style>
""", unsafe_allow_html=True)

# --- MAIN HERO SECTION ---
st.markdown('<div class="hero-title">In Silico Insecticide Discovery Pipeline</div>', unsafe_allow_html=True)
st.markdown('<div class="hero-subtitle">M.Sc. Dissertation Project • Genomics & Bioinformatics</div>', unsafe_allow_html=True)

st.divider()

# --- PROJECT OVERVIEW ---
st.write("### 📂 Project Modules")
st.write("Select a tool from the sidebar or the cards below to begin.")

# Create 3 columns for a "Card" layout
col1, col2 = st.columns(2)

with col1:
    st.markdown("""
    <div class="card">
        <div class="icon">🧬</div>
        <h3>Target Identification</h3>
        <p>Homology-based mining of transcriptomes to discover novel protein targets in pest species.</p>
        <p style="color: #64748B; font-size: 0.9em;">Features: BLAST Integration, ORF Finding, ESMFold Structure Prediction.</p>
    </div>
    """, unsafe_allow_html=True)
    # Button to guide user
    st.info("👈 Go to **01 Target Discovery** in the sidebar")

with col2:
    st.markdown("""
    <div class="card">
        <div class="icon">💊</div>
        <h3>Lead Molecule Screening</h3>
        <p>Chemoinformatics pipeline to find insect-specific ligands using ChEMBL data.</p>
        <p style="color: #64748B; font-size: 0.9em;">Features: Tanimoto Similarity, Lipinski Rules, Toxicity Filtering.</p>
    </div>
    """, unsafe_allow_html=True)
    st.info("👈 Go to **02 Insect Drug Search** in the sidebar")

# --- FOOTER ---
st.markdown("---")
st.markdown("""
**Student Name:** [Your Name]  
**Supervisor:** [Mentor Name]  
*Department of Life Sciences, [University Name]*
""")
