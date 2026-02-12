import streamlit as st
import streamlit.components.v1 as components
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import Entrez
import py3Dmol
import requests
from io import StringIO

# --- CONFIGURATION ---
st.set_page_config(layout="wide", page_title="In Silico Gene Discovery")
Entrez.email = "student.dissertation@university.edu"  # Generic academic email to satisfy NCBI

# --- CUSTOM CSS ---
st.markdown("""
<style>
    .stButton>button {
        width: 100%;
        background-color: #FF4B4B;
        color: white;
        font-weight: bold;
        border-radius: 8px;
        height: 3em;
    }
    .success-box {
        padding: 1rem;
        border-radius: 0.5rem;
        background-color: #d4edda;
        color: #155724;
        border: 1px solid #c3e6cb;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

# --- HELPER FUNCTIONS ---

def get_esmfold_structure(sequence):
    """Fetches 3D structure from ESMFold API."""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        # verify=False prevents SSL errors on college networks
        response = requests.post(url, data=sequence, verify=False) 
        if response.status_code == 200:
            return response.text
        return None
    except:
        return None

def render_mol_safe(pdb_str):
    """Renders molecule using pure HTML/JS (No fragile dependencies)."""
    view = py3Dmol.view(width=500, height=400)
    view.addModel(pdb_str, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    view_html = view._make_html()
    components.html(view_html, height=400, width=500)

def find_orf(sequence):
    """Finds the longest ORF starting with 'M'."""
    if isinstance(sequence, str):
        seq_obj = Seq(sequence)
    else:
        seq_obj = sequence
        
    frames = []
    # Translate all 3 forward frames
    for frame in range(3):
        try:
            remainder = len(seq_obj[frame:]) % 3
            if remainder > 0:
                trans = seq_obj[frame:-remainder].translate()
            else:
                trans = seq_obj[frame:].translate()
            frames.append(str(trans))
        except: continue
        
    longest_protein = ""
    for frame_seq in frames:
        parts = frame_seq.split("*")
        for part in parts:
            if "M" in part:
                # Find the first Methionine
                m_index = part.find("M")
                protein = part[m_index:]
                if len(protein) > len(longest_protein):
                    longest_protein = protein
    
    # If no ORF found (common in short fragments), return the longest translation frame
    if not longest_protein and frames:
        longest_protein = max(frames, key=len)
        
    return longest_protein

# --- MAIN APP LAYOUT ---

st.title("🧬 Homology-Based Gene Discovery")
st.markdown("### Transcriptome Mining & Structure Prediction")

col1, col2 = st.columns([1, 2])

with col1:
    st.info("💡 **Workflow:**\n1. Enter Sequence.\n2. Define Organism.\n3. BLAST finds gene fragments.\n4. Tool assembles & folds protein.")
    
    # INPUTS
    query_seq = st.text_area("1. Query Sequence (Protein or Nucleotide)", height=150, 
                             value=">sp|P07140|ACES_DROME Acetylcholinesterase\nMIRSSEVESGVGRVGVGVPGLVGVPGLVRLGSPQLRLRLLRLLRLRLRLRLRLRLRLRLRLRLRLGAGGAGGAGGAGGEGGAGGAGGAGGEGGAGGAGGLGPHPLPLLSLPLLPLLLLLGDGALSQEGRPEDGILEYVVRTVRGLRAPVCPASPNCQWPLTVTGVFGGRMQLVDTAALEQGVGQTPVESGALVRGIRQVLRLGSPWLPEKPTTPLGRDRVNCQWPLTVTGVFGGRMQLVDTAALEQGVGQTPVESGALVRGIRQVLRLGSPWLPEKPTTPLGRDRVNCQWPLTVTGVFGGRMQLVDTAALEQGVGQTPVESGALVRGIRQVLRLGSPWLPEKPTTPLGRDRVNCQWPL")
    
    target_organism = st.text_input("2. Target Organism", value="Helicoverpa armigera")
    
    program = st.selectbox("3. BLAST Program", ["tblastn (Protein Query)", "blastn (Nucleotide Query)"], index=0)
    
    st.caption("Note: 'tblastn' matches Protein input against DNA databases.")
    
    run_btn = st.button("🚀 Run Discovery Pipeline")

with col2:
    if run_btn:
        # --- INPUT VALIDATION ---
        if not query_seq.strip():
            st.error("⚠️ Please enter a query sequence.")
            st.stop()
            
        if not target_organism.strip():
            st.error("⚠️ Please enter a target organism (e.g., 'Musca domestica').")
            st.stop()

        status = st.empty()
        status.info(f"⏳ **Step 1:** BLAST searching against *{target_organism}*... (This uses remote NCBI servers, please wait)")
        
        try:
            # 1. PREPARE BLAST
            clean_query = query_seq
            if ">" in query_seq:
                clean_query = "".join([line for line in query_seq.splitlines() if not line.startswith(">")])
            
            # Determine program
            is_tblastn = "tblastn" in program
            prog_code = "tblastn" if is_tblastn else "blastn"
            
            organism_query = f"\"{target_organism}\"[Organism]"
            
            # 2. EXECUTE BLAST
            result_handle = NCBIWWW.qblast(
                prog_code, 
                "nt", 
                clean_query, 
                entrez_query=organism_query,
                expect=10.0,
                hitlist_size=5
            )
            
            # 3. PARSE RESULTS
            blast_record = NCBIXML.read(result_handle)
            
            if blast_record.alignments:
                top_alignment = blast_record.alignments[0]
                top_hsp = top_alignment.hsps[0]
                
                # Success Message
                st.markdown(f"""
                <div class="success-box">
                    <b>✅ Target Discovered!</b><br>
                    <b>Gene:</b> {top_alignment.title[:80]}...<br>
                    <b>E-value:</b> {top_hsp.expect}
                </div>
                """, unsafe_allow_html=True)
                
                # 4. PROCESS SEQUENCE
                st.markdown("### 🧬 Step 2: Sequence Assembly")
                
                # Get the sequence from the match
                # CLEANUP: Remove gaps (-) and numbers if any
                raw_seq = top_hsp.sbjct.replace("-", "").replace("\n", "")
                
                status.info("⏳ **Step 3:** Preparing Structure...")
                
                predicted_protein = ""
                
                # --- LOGIC FIX IS HERE ---
                if is_tblastn:
                    # tblastn returns PROTEIN. Use it directly.
                    st.info("ℹ️ TBLASTN used: Sequence is already translated.")
                    predicted_protein = raw_seq
                else:
                    # blastn returns DNA. We must translate it.
                    st.info("ℹ️ BLASTN used: Translating DNA to Protein...")
                    predicted_protein = find_orf(raw_seq)
                    # Fallback if ORF finder fails
                    if len(predicted_protein) < 10:
                        if len(raw_seq) % 3 != 0: raw_seq = raw_seq[:-(len(raw_seq)%3)]
                        predicted_protein = str(Seq(raw_seq).translate())

                st.write(f"**Predicted Protein Sequence ({len(predicted_protein)} AA):**")
                st.code(predicted_protein, language="text")
                
                # 5. FOLD STRUCTURE (ESMFold)
                st.markdown("### 🧊 Step 4: 3D Structure Prediction (ESMFold)")
                
                if len(predicted_protein) > 10:
                    pdb_data = get_esmfold_structure(predicted_protein)
                    
                    if pdb_data:
                        render_mol_safe(pdb_data)
                        st.download_button("📥 Download PDB File", pdb_data, "predicted_target.pdb")
                        status.success("✅ Pipeline Complete!")
                    else:
                        st.error("❌ Structure prediction failed (API might be busy or sequence too long).")
                else:
                    st.error("❌ Protein too short for folding.")

            else:
                status.error(f"❌ No matches found in '{target_organism}'.")
                st.warning("Suggestion: Try a broader organism query.")

        except Exception as e:
            st.error(f"❌ An error occurred: {e}")
