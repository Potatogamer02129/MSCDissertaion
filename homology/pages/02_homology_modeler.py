import streamlit as st
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio import Entrez
import pandas as pd
import requests
import py3Dmol
from stmol import showmol

# --- CONFIGURATION ---
st.set_page_config(layout="wide", page_title="In Silico Gene Discovery")
Entrez.email = "your.email@example.com"  # Always good practice to set this

# --- CUSTOM CSS FOR "AMAZING" LOOK ---
st.markdown("""
<style>
    .stButton>button {
        width: 100%;
        background-color: #FF4B4B;
        color: white;
        font-weight: bold;
    }
    .reportview-container {
        background: #f0f2f6;
    }
    h1 { color: #0E1117; }
    h3 { color: #262730; }
</style>
""", unsafe_allow_html=True)

# --- HELPER FUNCTIONS ---

def get_esmfold_structure(sequence):
    """Fetches 3D structure from ESMFold API (Faster than AlphaFold)."""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, verify=False) # verify=False for robustness in some networks
        if response.status_code == 200:
            return response.text
        else:
            return None
    except:
        return None

def find_orf(sequence):
    """Simple ORF finder: looks for longest protein between Start (M) and Stop (*)."""
    # Create a Seq object
    if isinstance(sequence, str):
        seq_obj = Seq(sequence)
    else:
        seq_obj = sequence
        
    # Translate all 3 forward frames
    frames = []
    for frame in range(3):
        trans = seq_obj[frame:].translate()
        frames.append(str(trans))
        
    # Find the longest continuous peptide starting with M
    longest_protein = ""
    for frame_seq in frames:
        parts = frame_seq.split("*")
        for part in parts:
            if "M" in part:
                # Find the first M
                m_index = part.find("M")
                protein = part[m_index:]
                if len(protein) > len(longest_protein):
                    longest_protein = protein
                    
    return longest_protein

# --- MAIN APP LAYOUT ---

st.title("🧬 Homology-Based Gene Discovery & Modeler")
st.markdown("### From Transcriptome to 3D Target Structure")

col1, col2 = st.columns([1, 2])

with col1:
    st.info("💡 **Mentor's Workflow:**\n1. Input Query Sequence.\n2. BLAST against Target Organism.\n3. Extract Nucleotide/Exons.\n4. Find ORF (Open Reading Frame).\n5. Predict 3D Structure.")
    
    # Input Section
    query_seq = st.text_area("Paste Query Sequence (Nucleotide or Protein)", height=150, placeholder=">Sequence_1\nATGC...")
    
    target_organism = st.text_input("Target Organism (Latin Name or TaxID)", placeholder="e.g., Helicoverpa armigera", value="Helicoverpa armigera")
    
    program = st.selectbox("BLAST Program", ["tblastn (Protein vs Trans. Genome)", "blastn (Nuc vs Nuc)"])
    
    st.caption("Note: 'tblastn' takes a protein and searches the target's DNA genome.")
    
    run_btn = st.button("🚀 Run Discovery Pipeline")

with col2:
    if run_btn and query_seq and target_organism:
        status = st.empty()
        
        # 1. BLAST SEARCH
        status.info(f"⏳ Blasting against **{target_organism}** database... This may take 1-2 mins.")
        
        try:
            # We use 'nr' (nucleotide collection) or 'wgs' if accessible, 
            # but for stability 'nt' filtered by organism is best for viva demonstration.
            if "tblastn" in program:
                blast_prog = "tblastn"
            else:
                blast_prog = "blastn"
                
            # Perform the BLAST
            result_handle = NCBIWWW.qblast(
                blast_prog, 
                "nt", 
                query_seq, 
                entrez_query=f"\"{target_organism}\"[Organism]"
            )
            
            blast_record = NCBIXML.read(result_handle)
            
            # 2. PARSE RESULTS
            if len(blast_record.alignments) > 0:
                top_alignment = blast_record.alignments[0]
                top_hsp = top_alignment.hsps[0] # Taking the best High Scoring Pair
                
                status.success("✅ Match Found!")
                
                st.subheader("🔍 Discovery Results")
                st.write(f"**Best Hit:** {top_alignment.title[:60]}...")
                st.write(f"**E-value:** {top_hsp.expect}")
                
                # Extract the hit sequence (The "Exon" or transcript part)
                hit_seq = top_hsp.sbjct
                
                # 3. ORF FINDING & TRANSLATION
                st.markdown("---")
                st.write("⚙️ **Processing Sequence & Finding ORF...**")
                
                # If tblastn, the hit is already translated in the alignment usually, 
                # but let's grab the DNA and re-translate to be safe/thorough for the "ORF" requirement.
                # For simplicity in this demo, we take the aligned subject sequence.
                
                # Clean gaps
                clean_nuc_seq = hit_seq.replace("-", "")
                final_protein = find_orf(clean_nuc_seq)
                
                if not final_protein:
                    # Fallback if ORF finder fails (e.g. if it's a short fragment)
                    final_protein = top_hsp.sbjct if blast_prog == "blastp" else clean_nuc_seq
                    st.warning("⚠️ No clear ORF starting with 'M' found. Using raw translation.")

                st.code(final_protein, language="text", line_numbers=True)
                st.caption(f"Predicted Protein Length: {len(final_protein)} AA")

                # 4. STRUCTURE PREDICTION (ESMFold)
                st.markdown("---")
                st.subheader("nPDB Structure Prediction (ESMFold)")
                
                if len(final_protein) > 400:
                    st.warning("⚠️ Protein is long. Rendering might be slow.")
                
                status.info("⏳ Generating 3D Structure...")
                pdb_string = get_esmfold_structure(final_protein)
                
                if pdb_string:
                    status.success("✅ Structure Generated!")
                    
                    # Render using stmol
                    view = py3Dmol.view(width=500, height=400)
                    view.addModel(pdb_string, 'pdb')
                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                    view.zoomTo()
                    showmol(view, height=400, width=500)
                    
                    # Download Button
                    st.download_button("📥 Download PDB File", pdb_string, "predicted_structure.pdb")
                else:
                    st.error("❌ Could not generate structure (API limit or sequence issue).")

            else:
                status.error("❌ No significant matches found in this organism.")
                
        except Exception as e:
            status.error(f"❌ An error occurred: {e}")
            st.error("Tip: Check your internet connection or try a different target organism.")
    
    elif run_btn:
        st.warning("⚠️ Please provide both a sequence and an organism name.")
