import streamlit as st
import streamlit.components.v1 as components
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import Entrez
import py3Dmol
import requests
import time

# --- 1. PAGE CONFIGURATION ---
st.set_page_config(
    layout="wide", 
    page_title="In Silico Target Discovery",
    page_icon="🧬"
)

# Academic email for NCBI (Best practice)
Entrez.email = "researcher.student@university.edu"

# --- 2. THEME-AGNOSTIC STYLING (CSS) ---
st.markdown("""
<style>
    /* 1. REMOVED hardcoded backgrounds (white/gray) for the main app.
       Streamlit will now automatically switch between Dark/Light mode 
       based on the user's browser settings.
    */
    
    /* 2. Success/Status Cards - We force these to be light colored with dark text
       so they pop out regardless of whether the app is in Dark or Light mode. */
    .success-card {
        padding: 1rem;
        border-radius: 8px;
        background-color: #DCFCE7; /* Light Green */
        border: 1px solid #86EFAC;
        color: #14532D; /* Dark Green Text (Always visible on light green) */
        margin-bottom: 1rem;
    }
    
    /* 3. Metric Containers - Use a subtle border instead of a solid white background 
       so it blends into both Dark and Light themes. */
    div[data-testid="metric-container"] {
        border: 1px solid rgba(128, 128, 128, 0.2);
        padding: 10px;
        border-radius: 8px;
        /* No background-color defined, so it takes the theme's color */
    }

    /* 4. Button Styling */
    .stButton>button {
        width: 100%;
        border-radius: 8px;
        height: 3em;
        font-weight: bold;
        /* Streamlit handles the button colors automatically */
    }
</style>
""", unsafe_allow_html=True)

# --- 3. ROBUST LOGIC FUNCTIONS (Cached for Speed) ---

@st.cache_data(show_spinner=False)
def run_blast_search(program, query_seq, organism):
    """
    Runs the BLAST search and returns the parsed record.
    Cached: If you search the same thing twice, it loads instantly.
    """
    # 1. Input Sanitization
    clean_query = query_seq
    if ">" in query_seq:
        clean_query = "".join([line for line in query_seq.splitlines() if not line.startswith(">")])
    
    prog_code = "tblastn" if "tblastn" in program else "blastn"
    organism_query = f"\"{organism}\"[Organism]"
    
    # 2. Execute Request
    result_handle = NCBIWWW.qblast(
        prog_code, 
        "nt", 
        clean_query, 
        entrez_query=organism_query,
        expect=10.0,
        hitlist_size=1  # We focus on the best hit for this tool
    )
    
    # 3. Parse and Return
    blast_record = NCBIXML.read(result_handle)
    return blast_record, prog_code

@st.cache_data(show_spinner=False)
def fold_protein_esm(sequence):
    """Fetches structure from ESMFold. Cached to save API limits."""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, verify=False)
        if response.status_code == 200:
            return response.text
        return None
    except:
        return None

def find_orf_robust(sequence):
    """
    Scans all 3 reading frames to find the longest reliable protein.
    """
    if isinstance(sequence, str):
        seq_obj = Seq(sequence)
    else:
        seq_obj = sequence
        
    frames = []
    # Translate 3 forward frames
    for frame in range(3):
        try:
            remainder = len(seq_obj[frame:]) % 3
            trim_seq = seq_obj[frame:-remainder] if remainder > 0 else seq_obj[frame:]
            frames.append(str(trim_seq.translate()))
        except: continue
        
    longest_protein = ""
    for frame_seq in frames:
        # Split by stop codons (*)
        parts = frame_seq.split("*")
        for part in parts:
            if "M" in part:
                # Start from Methionine
                protein = part[part.find("M"):]
                if len(protein) > len(longest_protein):
                    longest_protein = protein
    
    # Fallback: if no "Start-to-Stop" found, take longest generic fragment
    if not longest_protein and frames:
        longest_protein = max(frames, key=len)
        
    return longest_protein

def render_mol_advanced(pdb_str, style='cartoon', color='spectrum'):
    """Advanced molecule renderer with style options."""
    view = py3Dmol.view(width=700, height=500)
    view.addModel(pdb_str, 'pdb')
    
    if style == 'cartoon':
        view.setStyle({'cartoon': {'color': color}})
    elif style == 'stick':
        view.setStyle({'stick': {}})
    elif style == 'sphere':
        view.setStyle({'sphere': {'scale': 0.3}})
        
    view.zoomTo()
    view_html = view._make_html()
    components.html(view_html, height=500, width=700)

# --- 4. MAIN DASHBOARD UI ---

# Sidebar for Controls
with st.sidebar:
    st.header("⚙️ Search Parameters")
    
    with st.expander("📌 Search Instructions", expanded=True):
        st.write("""
        1. Paste a **known protein** (e.g., from Fruit Fly).
        2. Enter the **pest/insect** name you want to target.
        3. Click Run.
        """)
    
    target_organism = st.text_input("Target Organism", value="Helicoverpa armigera", help="Scientific name of the species")
    program = st.selectbox("Algorithm", ["tblastn (Protein → DNA)", "blastn (DNA → DNA)"])
    
    st.markdown("---")
    st.header("🧪 Visualization Settings")
    mol_style = st.selectbox("Structure Style", ["cartoon", "stick", "sphere"])
    
    st.markdown("---")
    st.caption("Dissertation Project Toolkit v1.0")

# Main Content Area
st.title("🧬 Homology-Based Target Identification")
st.markdown("""
This pipeline automates the **In Silico Cloning** process: identifying orthologous gene targets in non-model organisms and predicting their tertiary structure for druggability assessment.
""")

# Input Area
query_seq = st.text_area("Query Sequence (FASTA Format)", height=120, placeholder=">Sequence_ID\nMIRSSEVESGVGR...", help="Paste your reference protein or nucleotide sequence here.")

col1, col2, col3 = st.columns([1, 1, 2])
with col1:
    run_btn = st.button("🚀 Identify Target", help="Initiate remote BLAST search and structure prediction.")

# --- 5. EXECUTION LOGIC ---

if run_btn and query_seq:
    # A. Search Phase
    with st.status("Running Discovery Pipeline...", expanded=True) as status:
        
        st.write(f"1️⃣ Connecting to NCBI Remote Servers ({program})...")
        start_time = time.time()
        
        try:
            # Call Cached Function
            blast_record, used_prog = run_blast_search(program, query_seq, target_organism)
            
            if not blast_record.alignments:
                status.update(label="Search Failed", state="error", expanded=True)
                st.error(f"No significant homology found in *{target_organism}*.")
                st.stop()
                
            top_hit = blast_record.alignments[0]
            hsp = top_hit.hsps[0]
            
            # Calculate simple metrics
            ident = round((hsp.identities / hsp.align_length) * 100, 1)
            e_val = hsp.expect
            
            st.write(f"2️⃣ Match Found! (Identity: {ident}%, E-value: {e_val})")
            
            # B. Sequence Processing Phase
            st.write("3️⃣ Extracting & Assembling Sequence...")
            raw_seq = hsp.sbjct.replace("-", "").replace("\n", "")
            
            final_protein = ""
            if used_prog == "tblastn":
                final_protein = raw_seq # Already translated by BLAST
            else:
                final_protein = find_orf_robust(raw_seq)
                
            st.write(f"4️⃣ Folding Protein (Length: {len(final_protein)} AA)...")
            
            # C. Folding Phase
            if len(final_protein) > 10:
                pdb_data = fold_protein_esm(final_protein)
                if pdb_data:
                    status.update(label="Pipeline Completed Successfully!", state="complete", expanded=False)
                    
                    # --- 6. RESULTS DASHBOARD ---
                    
                    st.divider()
                    
                    # Summary Metrics Row
                    m1, m2, m3, m4 = st.columns(4)
                    m1.metric("Target Organism", target_organism)
                    m2.metric("Sequence Identity", f"{ident}%")
                    m3.metric("E-Value", f"{e_val:.1e}")
                    m4.metric("Protein Length", f"{len(final_protein)} AA")
                    
                    # Tabs for detailed view
                    tab1, tab2, tab3 = st.tabs(["🧊 3D Structure", "🧬 Sequences", "📋 Alignment Details"])
                    
                    with tab1:
                        st.markdown("#### Predicted Tertiary Structure (ESMFold)")
                        # Visual spacer
                        st.write("") 
                        render_mol_advanced(pdb_data, style=mol_style)
                        
                        # Download Button
                        st.download_button(
                            label="📥 Download PDB Structure",
                            data=pdb_data,
                            file_name=f"{target_organism.replace(' ', '_')}_target.pdb",
                            mime="text/plain"
                        )
                        
                    with tab2:
                        st.subheader("Recovered Target Sequence")
                        st.text_area("Predicted Protein", value=final_protein, height=150)
                        st.download_button("Copy/Download Fasta", f">Target_{target_organism}\n{final_protein}", "target_seq.fasta")
                        
                    with tab3:
                        st.subheader("BLAST Alignment Data")
                        st.write(f"**Definition:** {top_hit.title}")
                        st.write(f"**Accession:** {top_hit.accession}")
                        st.code(f"Query: {hsp.query[:80]}...\nMatch: {hsp.match[:80]}...\nSbjct: {hsp.sbjct[:80]}...", language="text")

                else:
                    status.update(label="Structure Prediction Failed", state="error")
                    st.error("ESM API Error: Sequence might be too complex for real-time folding.")
            else:
                status.update(label="Sequence Error", state="error")
                st.error("Recovered sequence is too short to be a valid protein.")
                
        except Exception as e:
            status.update(label="System Error", state="error")
            st.error(f"An unexpected error occurred: {str(e)}")
            st.info("Tip: Check your internet connection or the spelling of the organism.")
