import requests
from ..utils import fetch_json

def fetch_hpa_data(ensembl_id):
    """
    Fetches tissue expression data from Human Protein Atlas (HPA).
    Returns a list of dicts: {'tissue': str, 'value': float, 'unit': 'nTPM'}
    """
    if not ensembl_id:
        return []
        
    url = f"https://www.proteinatlas.org/{ensembl_id}.json"
    
    try:
        # HPA can be slow or flaky, simple get
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return []
            
        data = r.json()
        
        # We target "RNA tissue specific nTPM" as the high-confidence set
        # But for broader coverage (all detected), "RNA tissue specific nTPM" usually contains the enriched ones.
        # Let's check "RNA tissue specific nTPM" first.
        # Based on inspection, this is a dictionary: {"tissue": "value", ...}
        
        results = []
        
        # 1. Specific nTPM (Best quality)
        ntpm_dict = data.get("RNA tissue specific nTPM", {})
        if not ntpm_dict and "RNA tissue distribution" in data:
             # If specific is empty, maybe none are "Specific". 
             # We could fallback to "RNA tissue specific nTPM" being empty means no high expression?
             pass
             
        # Actually in the inspection script, ZAP70 (immune specific) had entries here.
        # If a protein is "Low tissue specificity" (expressed everywhere), this dict might be huge or empty?
        # Let's try to grab whatever is available in the nTPM dictionary.
        
        if isinstance(ntpm_dict, dict):
            for tissue, val_str in ntpm_dict.items():
                try:
                    val = float(val_str)
                    results.append({
                        "tissue": tissue,
                        "value": val,
                        "unit": "nTPM",
                        "source": "HPA"
                    })
                except:
                    continue
                    
        return results

    except Exception as e:
        print(f"Error HPA fetch for {ensembl_id}: {e}")
        return []
