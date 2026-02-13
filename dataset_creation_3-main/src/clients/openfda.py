import requests
from ..utils import fetch_json

OPENFDA_EVENT_URL = "https://api.fda.gov/drug/event.json"
OPENFDA_LABEL_URL = "https://api.fda.gov/drug/label.json"

def fetch_side_effects(drug_name, limit=50):
    """
    Fetches adverse events (side effects) for a given drug from OpenFDA.
    Returns a list of side effect terms.
    """
    if not drug_name:
        return []
        
    # Search for medicinalproduct matching the drug name
    # OpenFDA exact search is case sensitive? Try lowercase usually safe.
    query = f"patient.drug.medicinalproduct:{drug_name.lower()}"
    params = {
        "search": query,
        "count": "patient.reaction.reactionmeddrapt.exact",
        "limit": limit
    }
    
    try:
        data = fetch_json(OPENFDA_EVENT_URL, params=params)
        if data and "results" in data:
            # Extract terms from the count results
            return [item["term"] for item in data["results"]]
    except Exception as e:
        print(f"Error fetching side effects for {drug_name}: {e}")
    
    return []

def fetch_drug_label(drug_name):
    """
    Fetches the FDA Drug Label (package insert) text for a given drug.
    Returns a dict with key sections or None.
    """
    if not drug_name:
        return None

    # Search against the openfda.brand_name or generic_name
    # Using 'openfda.brand_name' or 'openfda.generic_name' is common.
    # We'll use a broad search on the indexed fields.
    query = f'openfda.brand_name:"{drug_name}" OR openfda.generic_name:"{drug_name}"'
    
    params = {
        "search": query,
        "limit": 1
    }
    
    try:
        data = fetch_json(OPENFDA_LABEL_URL, params=params)
        if data and "results" in data and len(data["results"]) > 0:
            res = data["results"][0]
            
            # Helper to join list of strings (OpenFDA returns lists for text blocks)
            def get_text(key):
                val = res.get(key)
                if isinstance(val, list):
                    return " ".join(val)
                return val

            return {
                # Rx vs OTC field mapping
                "mechanism_of_action": get_text("mechanism_of_action") or get_text("purpose") or get_text("active_ingredient"),
                "drug_interactions": get_text("drug_interactions"), # OTC usually puts this in 'warnings'
                "contraindications": get_text("contraindications") or get_text("do_not_use"),
                "boxed_warning": get_text("boxed_warning") or get_text("warnings"),
                "clinical_pharmacology": get_text("clinical_pharmacology") or get_text("indications_and_usage")
            }
            
    except Exception as e:
        print(f"Error fetching label for {drug_name}: {e}")
    
    return None
