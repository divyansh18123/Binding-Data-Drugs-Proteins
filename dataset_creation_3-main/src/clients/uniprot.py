from ..config import UNIPROT_BASE_URL, UNIPROT_JSON_URL
from ..utils import fetch_json
import urllib.parse

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

def fetch_uniprot_data(uniprot_id):
    url = UNIPROT_BASE_URL.format(uniprot_id)
    data = fetch_json(url)
    if data:
        return data
    # fallback: try explicit .json URL
    return fetch_json(UNIPROT_JSON_URL.format(uniprot_id))

def map_ensembl_to_uniprot(ensembl_ids):
    """
    Maps a list of Ensembl Protein IDs (ENSP...) to UniProt Accession IDs.
    Returns a dictionary {ensembl_id: uniprot_id}.
    """
    if not ensembl_ids:
        return {}
        
    # Clean IDs (strip 9606. prefix if present)
    clean_map = {}
    for eid in ensembl_ids:
        clean = eid.split(".")[-1] if "." in eid else eid
        clean_map[clean] = eid
        
    # Chunking to avoid URL length limits
    chunk_size = 10
    clean_ids = list(clean_map.keys())
    results = {}
    
    for i in range(0, len(clean_ids), chunk_size):
        chunk = clean_ids[i:i+chunk_size]
        # Query format: (xref:ensembl-ENSP000001) OR (xref:ensembl-ENSP000002)...
        query_parts = [f"(xref:ensembl-{cid})" for cid in chunk]
        query = " OR ".join(query_parts)
        
        params = {
            "query": query,
            "fields": "accession,xref_ensembl",
            "format": "json",
            "size": len(chunk) * 2 # Allow for some isoforms
        }
        
        data = fetch_json(UNIPROT_SEARCH_URL, params=params)
        
        if data and "results" in data:
            for item in data["results"]:
                acc = item.get("primaryAccession")
                # Find which input ID matches
                for ref in item.get("uniProtKBCrossReferences", []):
                    if ref.get("database") == "Ensembl":
                        # Check properties for ProteinId
                        # Format: id is Transcript (ENST), properties has ProteinId (ENSP)
                        candidate_ids = [ref.get("id")] # Add transcript just in case
                        
                        for prop in ref.get("properties", []):
                            if prop.get("key") == "ProteinId":
                                candidate_ids.append(prop.get("value"))
                        
                        for rid in candidate_ids:
                            # Strip version if present for matching (e.g. ENSP...5 -> ENSP...)
                            clean_rid = rid.split(".")[0]
                            
                            # Also check the raw rid just in case clean_map has it
                            if clean_rid in clean_map:
                                original_id = clean_map[clean_rid]
                                if original_id not in results:
                                    results[original_id] = acc
                                
    return results
