import streamlit as st

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="In Silico Pest Management Suite",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- THEME-AGNOSTIC CSS ---
st.markdown("""
<style>
    /* Gradient text for the main title */
    .main-title {
        background: -webkit-linear-gradient(45deg, #FF4B4B, #FF8F8F);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 3em;
        font-weight: bold;
        margin-bottom: 0px;
    }
    
    /* Subtitle styling */
    .subtitle {
        font-size: 1.2em;
        color: gray;
        margin-bottom: 2em;
    }
    
    /* Feature Card Styling - Adapts to Dark/Light mode */
    div[data-testid="column"] {
        border: 1px solid rgba(128, 128, 128, 0.2);
        border-radius: 10px;
        padding: 20px;
        transition: transform 0.2s;
    }
    div[data-testid="column"]:hover {
        transform: translateY(-5px);
        border: 1px solid #FF4B4B;
    }
</style>
""", unsafe_allow_html=True)

# --- SIDEBAR INFORMATION ---
with st.sidebar:
    st.image("https://cdn-icons-png.flaticon.com/512/2043/2043054.png", width=100) # Generic DNA/Bio icon
    st.header("MSc Dissertation Project")
    st.markdown("""
    **Developer:** Pratik Pareek  
    **Domain:** Bioinformatics & Drug Discovery  
    
    This platform integrates chemoinformatics and structural biology to identify novel insecticidal targets and predict hepatotoxicity.
    """)
    st.divider()
    st.caption("© 2026 | Built with Streamlit")

# --- MAIN PAGE CONTENT ---
st.markdown('<p class="main-title">In Silico Target Discovery Suite</p>', unsafe_allow_html=True)
st.markdown('<p class="subtitle">An integrated pipeline for novel insecticidal protein targeting and structural prediction.</p>', unsafe_allow_html=True)

st.write("Welcome to the project dashboard. Please select a module below or use the sidebar to navigate.")
st.write("") # Spacer

# --- TOOL CARDS ---
col1, col2, col3 = st.columns([1, 1, 0.2]) # Added an empty 3rd column for spacing

with col1:
    st.subheader("🧪 Module 1: Target Search")
    st.write("""
    **ChEMBL Database Mining** Search for the nearest possible insect protein targets using Tanimoto similarity and RDKit. Designed to overcome data scarcity in non-model organisms.
    """)
    st.write("") # Spacer
    # Streamlit's native page linker
    st.page_link("pages/01_ChEMBL_Target_Search.py", label="Open Target Search", icon="🔍")

with col2:
    st.subheader("🧬 Module 2: Homology Modeler")
    st.write("""
    **Transcriptome to Structure Pipeline** Perform cross-species gene discovery via NCBI BLAST, extract reading frames, and predict 3D tertiary structures in real-time using ESMFold.
    """)
    st.write("") # Spacer
    st.page_link("pages/02_Homology_Modeler.py", label="Open Homology Modeler", icon="🧊")

# --- DISSERTATION ABSTRACT/INFO ---
st.divider()
st.subheader("📌 Project Overview")
st.markdown("""
* **The Challenge:** Traditional pest control relies on broad-spectrum chemicals. Developing targeted solutions requires identifying species-specific proteins, but sequence data for agricultural pests is often fragmented.
* **The Solution:** This suite provides a two-pronged approach. Module 1 utilizes chemoinformatics to map known ligands to insect variants. Module 2 bridges the gap between raw transcriptomic data and actionable 3D structures, enabling rapid *in silico* cloning and structure-based drug design.
* **Future Work:** Integration of a predictive hepatotoxicity model for safety screening of proposed chemical agents.
""")
