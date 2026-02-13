from uniprot import fetch_uniprot_data
import quickgo
import uniprot
from tqdm import tqdm
import requests
import pandas as pd

import asyncio
import aiohttp

HEADERS_JSON = {"Accept": "application/json"}
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb/{}"
UNIPROT_JSON_URL  = "https://rest.uniprot.org/uniprotkb/{}.json"
QUICKGO_TERM_URL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}"
QUICKGO_ANCESTOR_URL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}/ancestors"
QUICKGO_MULTI_TERM_URL = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{}"
results_for_get_data = []
results_for_get_data2 = []
results_for_get_data3 = []


async def get_data(QUICKGO_TERM_URL, go_ids, results_for_get_data):
    async with aiohttp.ClientSession() as session:
        tasks = []
        for go_id in go_ids:
            tasks.append(asyncio.create_task(session.get(QUICKGO_TERM_URL.format(go_id), ssl=False)))
        responses = await asyncio.gather(*tasks)
        for response in responses:
            results_for_get_data.append(await response.json())
        

async def get_data2(QUICKGO_ANCESTOR_URL, go_ids_who_have_term):
    async with aiohttp.ClientSession() as session:
        tasks = []
        for go_id in go_ids_who_have_term:
            tasks.append(asyncio.create_task(session.get(QUICKGO_ANCESTOR_URL.format(go_id), ssl=False)))
        responses = await asyncio.gather(*tasks)
        for response in responses:
            results_for_get_data2.append(await response.json())





def fetch_json(url, headers=HEADERS_JSON, params=None, timeout=10):
    try:
        r = requests.get(url, headers=headers, params=params, timeout=timeout)
        r.raise_for_status()
        return r.json()
    except Exception:
        return None


# def fetch_go_parents_recursive_step(go_id):
#     # Implementation of parent fetching that matches main.py logic (single level or using ancestors)
#     # main.py actually fetched ancestors, then checked children of ancestors for relation.
#     # Here we can simplify or direct import that logic. 
#     # For now, I'll basically replicate main.py's specific logic for edge finding.
    
#     edges = []

#     data = fetch_json(QUICKGO_ANCESTOR_URL.format(go_id)) #*****
#     if data and data.get("results"):
#         ancestors = data["results"][0].get("ancestors", [])
#     else:
#         ancestors = []


#     if not ancestors:
#         return edges
    
#     # We need to find which ancestor is a direct parent.
#     # This part of main.py was a bit heavy (multi call).
#     # We can implement it in quickgo.py or here.
#     # Let's rely on quickgo helper.
    
#     # Actually logic in main.py:
#     # 1. get ancestors
#     # 2. get details of ancestors
#     # 3. iterate ancestor children to see if CHILD == original GO_ID
    
#     ancestors = [a for a in ancestors if a != go_id]
#     if not ancestors:
#         return edges

#     # Fetch details for all ancestors
#     # terms_data = quickgo.fetch_multiple_terms(ancestors)

#     if not ancestors:
#         terms_data = []
#     ids_str = ",".join(ancestors)
#     data = fetch_json(QUICKGO_MULTI_TERM_URL.format(ids_str)) #*****
#     if data:
#         terms_data = data.get("results", [])
#     else:
#         terms_data = []


#     for term_data in terms_data:
#         parent_id = term_data["id"]
#         children = term_data.get("children", [])
#         for child in children:
#             if child.get("id") == go_id:
#                 edges.append({
#                     "child_go_id": go_id,
#                     "parent_go_id": parent_id,
#                     "relation": child.get("relation"),
#                 })
#     return edges


# def fetch_uniprot_data(uniprot_id):
#     return fetch_json(UNIPROT_BASE_URL.format(uniprot_id))

# davis_df = pd.read_csv("davis.txt", sep=r"\t" ,engine = "python")
# kiba_df = pd.read_csv("KIBA_Unique_Proteins.xlsx - Sheet1.csv", engine = 'python')
# uniprot_ids = kiba_df["ProteinID"].unique()

# uniprot_ids = davis_df['uniprot_id'].unique()
# # uniprot_ids = ["P00519", "Q9BQB6"]
# df = pd.read_csv("metz.txt", engine = 'python', sep = '\t')
# Add the uniprot ids
uniprot_ids = df['uniprot_id'].unique()

go_term_cache = {}
data_store = {
            "ontology_edges": [],
            "ontology_terms": []
        }

for uid in tqdm(uniprot_ids, desc="Processing Proteins"):
    # print(uid)
    go_id_who_have_term = []
    data = uniprot.fetch_uniprot_data(uid)

    url = UNIPROT_BASE_URL.format(uid)
    data = fetch_json(url)
    if not data:
        data = fetch_json(UNIPROT_JSON_URL.format(uid))
        if not data:
            continue


    if not data:
        continue

    go_ids = set()
    pid = data.get("primaryAccession")

    for ref in data.get("uniProtKBCrossReferences", []):
        if ref["database"] == "GO":
            gid = ref["id"]
            go_ids.add(gid)
    


    asyncio.run(get_data(QUICKGO_TERM_URL, go_ids, results_for_get_data))


    for result in results_for_get_data:
        # term = quickgo.fetch_term(go_id)
        
        data = result
        if data and data.get("results"):
            term = data["results"][0]
        else:
            term = None
        
        go_id = term.get("id")
        
        if term:
            go_term_cache[go_id] = {
                "go_id": go_id,
                "name": term.get("name"),
                "aspect": term.get("aspect"),
                "definition": term.get("definition", {}).get("text"),
            }
            go_id_who_have_term.append(go_id)

            # parents_data = fetch_go_parents_recursive_step(go_id)
            # data_store["ontology_edges"].extend(parents_data)
    # asyncio.run(go_id_who_have_term)
    # asyncio.run(get_data2(QUICKGO_ANCESTOR_URL, go_id_who_have_term))
    
    # for result in results_for_get_data2:
    #     edges = []
    #     data = result #*****

    #     if data and data.get("results"):
    #         go_id = data["results"][0].get("id")
    #     else:
    #         go_id = None

    #     if data and data.get("results"):
    #         ancestors = data["results"][0].get("ancestors", [])
    #     else:
    #         ancestors = []


    #     if not ancestors:
    #         # data_store["ontology_edges"].extend([])
    #         continue
        
    #     ancestors = [a for a in ancestors if a != go_id]
    #     if not ancestors:
    #         # data_store["ontology_edges"].extend([])
    #         continue

    #     # Fetch details for all ancestors
    #     # terms_data = quickgo.fetch_multiple_terms(ancestors)
    #     print(result)
    #     print(ancestors)


    #     asyncio.run(get_data(QUICKGO_TERM_URL, ancestors, results_for_get_data3))

    #     for result in results_for_get_data3:

    #         if not ancestors:
    #             terms_data = []

    #         # ids_str = ",".join(ancestors)

    #         data = result #*****
    #         if data:
    #             terms_data = data.get("results", [])
    #         else:
    #             terms_data = []

    #         for term_data in terms_data:
    #             parent_id = term_data["id"]
    #             children = term_data.get("children", [])
    #             for child in children:
    #                 if child.get("id") == go_id:
    #                     edges.append({
    #                         "child_go_id": go_id,
    #                         "parent_go_id": parent_id,
    #                         "relation": child.get("relation"),
    #                     })
                    
    #     results_for_get_data3 = []
        
    #     data_store["ontology_edges"].extend(edges)


       
    results_for_get_data = []
    # results_for_get_data2 = []
    pd.DataFrame(go_term_cache.values()).to_csv("ontology_terms.csv", index = False)

# File 1 above

# # import requests
# # from .config import HEADERS_JSON

# # def fetch_json(url, headers=HEADERS_JSON, params=None, timeout=10):
# #     try:
# #         r = requests.get(url, headers=headers, params=params, timeout=timeout)
# #         r.raise_for_status()
# #         return r.json()
# #     except Exception:
# #         return None

# # def fetch_text(url, timeout=10):
# #     try:
# #         r = requests.get(url, timeout=timeout)
# #         if r.status_code == 200:
# #             return r.text
# #         return None
# #     except Exception:
# #         return None


# # import pandas as pd
# # pd.DataFrame(go_term_cache.values()).to_csv("ontology_terms.csv", index = False)
# # pd.DataFrame(data_store["ontology_edges"]).to_csv("ontology_edges.csv", index = False)

# ## For file 2

go_term_cache = {}
data_store = {
            "ontology_edges": [],
            "ontology_terms": []
        }

        
terms_combined_df = pd.read_csv("ontology_terms.csv", engine = "python")
# terms_combined_df = pd.read_csv("temp.csv", engine = "python")
go_ids = terms_combined_df["go_id"].unique()


asyncio.run(get_data2(QUICKGO_ANCESTOR_URL, go_ids))

for result in tqdm(results_for_get_data2, desc = "Go id processing: "):
    edges = []
    data = result #*****

    if data and data.get("results"):
        go_id = data["results"][0].get("id")
    else:
        go_id = None

    if data and data.get("results"):
        ancestors = data["results"][0].get("ancestors", [])
    else:
        ancestors = []


    if not ancestors:
        # data_store["ontology_edges"].extend([])
        continue
    
    ancestors = [a for a in ancestors if a != go_id]
    if not ancestors:
        # data_store["ontology_edges"].extend([])
        continue

    # Fetch details for all ancestors
    # terms_data = quickgo.fetch_multiple_terms(ancestors)


    try:
        asyncio.run(get_data(QUICKGO_TERM_URL, ancestors, results_for_get_data3))
    except:
        print("error", uid)

    for result in results_for_get_data3:

        if not ancestors:
            terms_data = []

        # ids_str = ",".join(ancestors)

        data = result #*****
        if data:
            terms_data = data.get("results", [])
        else:
            terms_data = []

        for term_data in terms_data:
            parent_id = term_data["id"]
            children = term_data.get("children", [])
            for child in children:
                if child.get("id") == go_id:
                    edges.append({
                        "child_go_id": go_id,
                        "parent_go_id": parent_id,
                        "relation": child.get("relation"),
                    })
                
    results_for_get_data3 = []
        
    data_store["ontology_edges"].extend(edges)
    pd.DataFrame(data_store["ontology_edges"]).to_csv("ontology_edges.csv", index = False)
