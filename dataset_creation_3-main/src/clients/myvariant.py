import requests
from ..utils import fetch_json

MYVARIANT_URL = "https://myvariant.info/v1/variant"

def fetch_variants_metadata(rsids):
    """
    Fetches detailed metadata for a list of rsIDs from MyVariant.info.
    Uses POST for batch retrieval.
    """
    if not rsids:
        return {}
    
    # MyVariant allows batch queries via POST (max 1000 usually)
    BATCH_SIZE = 1000
    mapping = {}
    
    unique_rsids = list(set(rsids))
    
    for i in range(0, len(unique_rsids), BATCH_SIZE):
        batch = unique_rsids[i:i + BATCH_SIZE]
        payload = {
            "ids": ",".join(batch),
            "fields": "clinvar.rcv.clinical_significance,clinvar.clinical_significance,snpeff.ann.effect,snpeff.ann.hgvs_p,dbsnp.rsid"
        }
        
        headers = { "Content-Type" : "application/x-www-form-urlencoded" }
        
        try:
            resp = requests.post(MYVARIANT_URL, data=payload, headers=headers)
            resp.raise_for_status()
            results = resp.json()
            
            for r in results:
                rsid = r.get("query")
                if not rsid: continue
                
                # Extract Clinical Significance (robustly)
                cs = r.get("clinvar", {}).get("clinical_significance")
                if not cs:
                    # Try RCV list
                    rcv = r.get("clinvar", {}).get("rcv", [])
                    if isinstance(rcv, dict): rcv = [rcv]
                    sig_set = set()
                    for entry in rcv:
                        sig = entry.get("clinical_significance")
                        if sig: sig_set.add(sig)
                    if sig_set:
                        cs = "; ".join(sorted(sig_set))
                
                mapping[rsid] = {
                    "clinical_significance": cs,
                    "consequence": _extract_consequence(r),
                    "protein_change": _extract_protein_change(r)
                }
        except Exception as e:
            print(f"Error fetching MyVariant batch: {e}")
            
    return mapping

def _extract_consequence(r):
    # snpeff.ann can be a list or dict
    ann = r.get("snpeff", {}).get("ann", [])
    if isinstance(ann, dict): ann = [ann]
    if ann:
        return ann[0].get("effect")
    return None

def _extract_protein_change(r):
    ann = r.get("snpeff", {}).get("ann", [])
    if isinstance(ann, dict): ann = [ann]
    if ann:
        return ann[0].get("hgvs_p")
    return None
