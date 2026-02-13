import requests
import time
from ..utils import fetch_json

STRING_API_URL = "https://string-db.org/api/json/network"

def fetch_interactions(protein_list, species=9606, min_score=400):
    """
    Fetch interactions for a list of proteins from STRING DB.
    
    Args:
        protein_list (list): List of protein identifiers (Gene names or Ensembl IDs).
        species (int): NCBI Taxon ID (default 9606 for Human).
        min_score (int): Minimum combined score (0-1000). Default 400 (Medium confidence).
        
    Returns:
        list: List of interaction dictionaries.
    """
    if not protein_list:
        return []
        
    # STRING API accepts 'identifiers' as newline separated string
    # Limit batch size to reasonable amount (e.g. 100)
    # For now, we assume this is called per protein or small batches.
    
    params = {
        "identifiers": "%0d".join(protein_list),
        "species": species,
        "caller_identity": "dataset_enrichment_project",
        "required_score": min_score
    }
    
    try:
        data = fetch_json(STRING_API_URL, params=params)
        return data if isinstance(data, list) else []
    except Exception as e:
        print(f"Error fetching STRING interactions: {e}")
        return []
