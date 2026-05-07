import pandas as pd
import requests
import sys
import os
import tqdm

# Add current directory to path to allow importing from src
sys.path.append(os.getcwd())

try:
    from src.clients import uniprot, pdb
except ImportError:
    print("Error: Could not import 'src.clients'. Make sure you are running this script from the root directory.")
    sys.exit(1)


def fetch_alphafold_entry(uniprot_id):
    """
    Fetches the AlphaFold DB entry ID and CIF URL for a given UniProt ID.
    Returns (entry_id, 'AlphaFoldDB', cif_url) or (None, None, None).
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        # Timeout is important as this is external
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data and isinstance(data, list):
                # Take the first entry
                entry = data[0]
                return entry.get('entryId'), "AlphaFoldDB", entry.get('cifUrl')
    except Exception as e:
        print(f"Error fetching AF for {uniprot_id}: {e}")
    return None, None, None

def get_best_structure(uniprot_id):
    """
    Fetches the best structure ID for a given UniProt ID.
    Priority:
    1. Experimental PDB structure (lowest resolution)
    2. AlphaFoldDB structure (fetched via API)
    
    Returns: (structure_id, structure_source, structure_url) or (None, None, None)
    """
    try:
        data = uniprot.fetch_uniprot_data(uniprot_id)
        if not data:
            return None, None, None
        
        # 1. Look for PDBs
        pdb_ids = set()
        has_alphafold_ref = False
        
        for ref in data.get("uniProtKBCrossReferences", []):
            if ref["database"] == "PDB":
                pdb_ids.add(ref["id"])
            elif ref["database"] == "AlphaFoldDB":
                has_alphafold_ref = True
        
        # If PDBs exist, verify resolution and pick best
        if pdb_ids:
            best_pdb = None
            best_resolution = float('inf')
            
            for pid in pdb_ids:
                pdb_data = pdb.fetch_pdb_entry(pid)
                if not pdb_data:
                    continue
                    
                # Get resolution
                res_list = pdb_data.get("rcsb_entry_info", {}).get("resolution_combined", [])
                resolution = res_list[0] if res_list else None
                
                if resolution is not None:
                    if resolution < best_resolution:
                        best_resolution = resolution
                        best_pdb = pid
                elif best_pdb is None:
                    # Provide a fallback if no resolution found but ID exists
                    best_pdb = pid
            
            if best_pdb:
                # Construct RCSB CIF URL
                cif_url = f"https://files.rcsb.org/download/{best_pdb}.cif"
                return best_pdb, "PDB", cif_url

        # 2. Fallback to AlphaFold if reference exists or try anyway
        if has_alphafold_ref:
            return fetch_alphafold_entry(uniprot_id)

        return None, None, None

    except Exception as e:
        print(f"Error fetching structure for {uniprot_id}: {e}")
        return None, None, None

def main():
    input_file = "davis_392.txt"
    # We will read the PREVIOUS output if it exists to save time
    existing_output_file = "davis_392_with_pdb.csv"
    output_file = "davis_392_with_pdb.csv"
    
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file, sep="\t") 
    
    # Load existing results if available
    existing_map = {} # uid -> (id, source, url)
    if os.path.exists(existing_output_file):
        print(f"Found existing {existing_output_file}, loading to cache...")
        try:
            existing_df = pd.read_csv(existing_output_file)
            
            # Use strict columns
            cols = existing_df.columns
            if 'structure_id' in cols and 'structure_source' in cols:
                # We optionally check for structure_url
                has_url = 'structure_url' in cols
                
                temp_cols = ['uniprot_id', 'structure_id', 'structure_source']
                if has_url: temp_cols.append('structure_url')
                
                temp = existing_df[temp_cols].drop_duplicates()
                
                for _, row in temp.iterrows():
                     uid = row['uniprot_id']
                     sid = row['structure_id']
                     source = row['structure_source']
                     url = row['structure_url'] if has_url else None
                     
                     if pd.notna(uid) and pd.notna(sid):
                         existing_map[uid] = (sid, source, url)
                         
            print(f"Loaded {len(existing_map)} existing mappings.")
        except Exception as e:
            print(f"Could not load existing file: {e}")

    unique_proteins = df['uniprot_id'].dropna().unique()
    print(f"Found {len(unique_proteins)} unique proteins.")
    
    structure_map = {} # uid -> (id, source, url)
    
    print("Fetching/Generating structure URLs...")
    for uid in tqdm.tqdm(unique_proteins):
        if pd.isna(uid): continue
        
        # Check cache/existing
        if uid in existing_map:
            cached = existing_map[uid]
            sid, source, url = cached
            
            # If we have a URL, great.
            if url and pd.notna(url):
                structure_map[uid] = cached
                continue
            
            # If we don't have a URL, optimize generation
            if source == 'PDB':
                new_url = f"https://files.rcsb.org/download/{sid}.cif"
                structure_map[uid] = (sid, source, new_url)
            elif source == 'AlphaFoldDB':
                # Call AF API to get versioned URL
                af_id, _, af_url = fetch_alphafold_entry(uid)
                if af_url:
                     structure_map[uid] = (af_id, source, af_url)
                else:
                     # Retry full logic if API fails
                     structure_map[uid] = get_best_structure(uid)
            else:
                structure_map[uid] = get_best_structure(uid)
        else:
            structure_map[uid] = get_best_structure(uid)
        
    # Map results
    print("Mapping IDs to DataFrame...")
    
    def get_id(uid):
        return structure_map.get(uid, (None, None, None))[0]
    
    def get_source(uid):
        return structure_map.get(uid, (None, None, None))[1]

    def get_url(uid):
        return structure_map.get(uid, (None, None, None))[2]

    df['structure_id'] = df['uniprot_id'].apply(get_id)
    df['structure_source'] = df['uniprot_id'].apply(get_source)
    df['structure_url'] = df['uniprot_id'].apply(get_url)
    df['pdb_id'] = df['structure_id'] 

    print(f"Saving to {output_file}...")
    df.to_csv(output_file, index=False)
    
    # Stats
    total = len(unique_proteins)
    found = df.drop_duplicates('uniprot_id')['structure_id'].notna().sum()
    with_url = df.drop_duplicates('uniprot_id')['structure_url'].notna().sum()
    pdb_count = df.drop_duplicates('uniprot_id')[df['structure_source'] == 'PDB']['structure_id'].count()
    af_count = df.drop_duplicates('uniprot_id')[df['structure_source'] == 'AlphaFoldDB']['structure_id'].count()
    
    print(f"Done.")
    print(f"Total Proteins: {total}")
    print(f"Mapped: {found} ({found/total:.1%})")
    print(f"URLs Generated: {with_url}")
    print(f"  - PDB (Experimental): {pdb_count}")
    print(f"  - AlphaFoldDB (Modeled): {af_count}")

if __name__ == "__main__":
    main()
