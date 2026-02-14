import pandas as pd
from chembl_webresource_client.new_client import new_client
import time

print("🧬 Initiating Automated ChEMBL Insecta Mining...")

# 1. Setup API Clients
target_api = new_client.target
activity_api = new_client.activity

# 2. Find all Insect Targets
print("🔍 Searching for all targets in taxonomy 'Insecta'...")
# Note: We filter for single proteins to ensure high-quality target data
insect_targets = target_api.filter(target_organism__icontains="Drosophila").filter(target_type="SINGLE PROTEIN")
# (Using Drosophila as a primary model to keep the API pull fast, 
#  but you can change to "Musca", "Anopheles", etc.)

target_dict = {}
for t in insect_targets:
    target_dict[t['target_chembl_id']] = {
        'Target_Name': t['pref_name'],
        'Target_UniProt': t.get('target_components', [{}])[0].get('accession', 'Unknown'),
        'Species': t['organism']
    }

print(f"✅ Found {len(target_dict)} high-quality insect protein targets.")

# 3. Pull Bioactive Ligands for these targets
print("⏳ Pulling active compounds (IC50/EC50) for these targets... (This may take a few minutes)")
dataset = []

for target_id, info in target_dict.items():
    # Fetch activities with strong binding affinity (IC50 < 10,000 nM)
    activities = activity_api.filter(target_chembl_id=target_id, standard_type__in=['IC50', 'EC50'], standard_value__lt=10000)
    
    count = 0
    for act in activities:
        # We only want entries that have a valid SMILES structure
        if act.get('canonical_smiles'):
            dataset.append({
                'Ligand_SMILES': act['canonical_smiles'],
                'Target_UniProt': info['Target_UniProt'],
                'Target_Name': info['Target_Name'],
                'Species': info['Species'],
                'Binding_Affinity_nM': act['standard_value'],
                'Metric': act['standard_type']
            })
            count += 1
            if count > 50:  # Cap at 50 ligands per target for dataset balance
                break
    
    print(f"   - Pulled {count} ligands for {info['Target_Name']}")
    time.sleep(1) # Be polite to the ChEMBL API

# 4. Save to CSV
df = pd.DataFrame(dataset)
df.drop_duplicates(subset=['Ligand_SMILES', 'Target_UniProt'], inplace=True)
df.to_csv("insect_ligand_targets_raw.csv", index=False)

print(f"\n🎉 Success! Dataset generated with {len(df)} curated insect ligands.")
print("Saved as 'insect_ligand_targets_raw.csv'. Your Streamlit app will now use this enriched data.")
