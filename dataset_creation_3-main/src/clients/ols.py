import requests
from ..utils import fetch_json
from urllib.parse import quote

OLS_API_URL = "https://www.ebi.ac.uk/ols4/api"

def fetch_term(ontology, term_id):
    """
    Fetches term details from OLS.
    """
    # Double encode the IRI if needed, but for OBO IDs (UBERON:0000001) usually replacing : with _ works for OLS
    # or passing the short_form.
    
    # UBERON:0000955 -> uberon/terms?iri=... or terms/UBERON:0000955
    # The OLS API prefers double encoded URIs or short forms.
    
    # Try searching for the term by short form or id
    url = f"{OLS_API_URL}/ontologies/{ontology}/terms"
    params = {"short_form": term_id}
    
    try:
        data = fetch_json(url, params=params)
        if data and "_embedded" in data and "terms" in data["_embedded"]:
            terms = data["_embedded"]["terms"]
            if terms:
                return terms[0] # Return best match
    except Exception as e:
        print(f"Error fetching OLS term {term_id}: {e}")
        
    return None

def fetch_uberon_term(uberon_id):
    term = fetch_term("uberon", uberon_id)
    if term:
        return {
            "id": uberon_id,
            "name": term.get("label"),
            "description": term.get("description", [""])[0] if term.get("description") else ""
        }
    return None

def search_uberon_by_name(name):
    """
    Searches UBERON for a tissue name. Returns (id, name, description, iri).
    """
    url = f"{OLS_API_URL}/search"
    params = {
        "q": name,
        "ontology": "uberon",
        "exact": "true", # Try exact first
        "queryFields": "label"
    }
    
    try:
        data = fetch_json(url, params=params)
        docs = data.get("response", {}).get("docs", [])
        
        # If no exact match, try loose
        if not docs:
            params["exact"] = "false"
            data = fetch_json(url, params=params)
            docs = data.get("response", {}).get("docs", [])
            
        if docs:
            top = docs[0]
            desc = top.get("description", [])
            if isinstance(desc, list) and desc: 
                desc_text = desc[0]
            else: 
                desc_text = ""
                
            return {
                "id": top.get("obo_id"),
                "name": top.get("label"),
                "description": desc_text,
                "iri": top.get("iri")
            }
    except Exception as e:
        print(f"Error searching OLS for {name}: {e}")
        
    return None

def fetch_term_parents(term_iri):
    """
    Fetches the direct parents of a term via its IRI.
    Returns list of {id, name, description}.
    """
    # Requires double encoding of IRI?
    # Usually we can look up the term by IRI then follow 'hierarchicalParents' link
    # But we can also use the /ontologies/uberon/terms endpoint with iri param
    
    url = f"{OLS_API_URL}/ontologies/uberon/terms"
    params = {"iri": term_iri}
    
    try:
        data = fetch_json(url, params=params)
        terms = data.get("_embedded", {}).get("terms", [])
        if not terms: return []
        
        main_term = terms[0]
        links = main_term.get("_links", {})
        
        if "hierarchicalParents" in links:
            parent_url = links["hierarchicalParents"]["href"]
            p_data = fetch_json(parent_url)
            parents = p_data.get("_embedded", {}).get("terms", [])
            
            results = []
            for p in parents:
                desc = p.get("description", [])
                desc_text = desc[0] if desc else ""
                results.append({
                    "id": p.get("obo_id"),
                    "name": p.get("label"),
                    "description": desc_text
                })
            return results
            
    except Exception as e:
        print(f"Error fetching OLS parents for {term_iri}: {e}")
        
    return []
