import requests
from ..utils import fetch_json

ENSEMBL_REST_URL = "https://rest.ensembl.org"

def fetch_variants(protein_id):
    """
    Fetches genetic variants for a protein (UniProt ID) from Ensembl.
    Note: Ensembl often maps UniProt ID to Gene/Transcript. 
    Here we try to map UniProt -> Genomic Coordinates or use overlap endpoint.
    
    A more direct way for protein variants is using the 'overlap' endpoint with the translation ID,
    but mapping UniProt -> Ensembl Protein ID (ENSP) is required first.
    
    For simplicity in this demo, we will use the UniProt features which often contain variants,
    OR we can try to look up the gene symbol in Ensembl.
    
    However, UniProt's own 'features' endpoint (already used) contains 'Natural Variant' data.
    We will strictly use this client if we want *genomic* context or dbSNP IDs that UniProt might miss.
    
    Let's implement a VEP (Variant Effect Predictor) style lookup or Phenotype lookup if possible.
    Actually, simpler: Map UniProt -> Gene -> Phenotypes.
    """
    # 1. Map UniProt ID to Ensembl Gene ID
    xref_url = f"{ENSEMBL_REST_URL}/xrefs/symbol/homo_sapiens/{protein_id}?external_db=Uniprot/SWISSPROT"
    # This is tricky because protein_id is not a symbol.
    
    # Better: Use the 'lookup' endpoint if we had the Ensembl ID.
    # Since we don't, we might rely on UniProt's cross-refs to Ensembl.
    
    # Alternative: Fetch phenotypes for the gene associated with the protein.
    # We will need the Gene Name for this.
    pass

def fetch_variants_by_gene(gene_name):
    """
    Fetches variants/phenotypes linked to a gene symbol.
    """
    if not gene_name:
        return []
        
    ext_url = f"{ENSEMBL_REST_URL}/phenotype/gene/homo_sapiens/{gene_name}?include_associated=1"
    headers = { "Content-Type" : "application/json"}
    
    try:
        data = fetch_json(ext_url, headers=headers)
        if isinstance(data, list):
            return data
    except Exception as e:
        print(f"Error fetching variants for {gene_name}: {e}")
        
    return []
