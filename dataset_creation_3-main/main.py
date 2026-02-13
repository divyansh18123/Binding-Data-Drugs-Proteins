import os
from src.processors.protein import ProteinProcessor
from src.processors.drug import DrugProcessor
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from src.clients import chembl, openfda, mychem
from tqdm import tqdm
import gc
from src.csv_writer import safe_to_csv

# threading below
# DATA_DIR = "data"
# MAX_DRUG_THREADS = 6
# MAX_PROTEIN_THREADS = 4

os.makedirs("data", exist_ok=True)

# def process_drug_worker(cid: str, write_header: bool):
#     drugProc = DrugProcessor()
#     data_store = drugProc.process_single_drug(cid)

#     if not data_store:
#         return []

#     for name, rows in data_store.items():
#         df = pd.DataFrame(rows)
#         safe_to_csv(
#             df,
#             f"{DATA_DIR}/{name}.csv",
#             header=write_header,
#             mode="a"
#         )

#     discovered_enzymes = drugProc.get_discovered_enzyme_ids()

#     del data_store
#     gc.collect()

#     return discovered_enzymes


# def process_protein_worker(uid: str, write_header: bool):
#     proteinProc = ProteinProcessor()
#     result = proteinProc.process_single_protein(uid, fetch_interactions=True)

#     if not result:
#         return

#     data_store, metadata = result

#     for name, rows in data_store.items():
#         df = pd.DataFrame(rows)
#         safe_to_csv(
#             df,
#             f"{DATA_DIR}/{name}.csv",
#             header=write_header,
#             mode="a"
#         )

#     for name, df in metadata.items():
#         safe_to_csv(
#             df,
#             f"{DATA_DIR}/{name}.csv",
#             header=write_header,
#             mode="a"
#         )

#     del data_store, metadata
#     gc.collect()




def main():
    # Ensure data directory exists
    
    '''    
    # 1. Process Drugs FIRST (to discover enzymes)
    # Imatinib (CHEMBL941), Warfarin (CHEMBL1459), Aspirin (CHEMBL25)
    # target_drugs = ["CHEMBL941"]
    # davis_df = pd.read_csv("davis.txt", sep=r"\t", engine="python")

    # kiba_df = pd.read_csv("KIBA_Unique_Compounds.xlsx - Sheet1.csv", engine = 'python')
    # drug_df = pd.read_csv("data_kiba/drug.csv", engine = 'python')
    # # print(protein_df["chembl_id"])
    # print(list(set(kiba_df["CHEMBLID"]) - set(drug_df["chembl_id"])))

    

    # target_drugs = kiba_df["CHEMBLID"].unique()
    
    # davis_df = pd.read_csv("davis.txt", sep=r"\t", engine="python")
    # target_drugs = davis_df["chembl_id"].unique()

    # target_drugs = target_drugs.unique()
    # target_drugs = ["CHEMBL396377"]

    # failed for mychem, list error kiba
    # target_drugs = ["CHEMBL7929", "CHEMBL279405", "CHEMBL265502", "CHEMBL260135", "CHEMBL1971694", "CHEMBL60254", "CHEMBL75232", "CHEMBL493937", "CHEMBL31", "CHEMBL1971289", "CHEMBL396377", "CHEMBL163269", "CHEMBL444864"]
    
    # print("=== Processing Drugs ===")
    # drug_proc = DrugProcessor()
    # drug_proc.process_drugs(target_drugs)
    
    # # Get discovered enzymes (UniProt IDs)
    # discovered_enzymes = drug_proc.get_discovered_enzyme_ids()
    # print(discovered_enzymes)
    # print(f"Discovered {len(discovered_enzymes)} enzymes from drug metabolism.")
    

    
    # 2. Process Proteins (Targets + Enzymes)
    # ABL1 (P00519) -> Target of Imatinib
    # VKORC1 (Q9BQB6) -> Target of Warfarin
    # PTGS1 (P23219) -> Target of Aspirin
    # base_proteins = ["P00519"] 
    # davis_df = pd.read_csv("davis.txt", sep=r"\t", engine="python")

    # kiba_df = pd.read_csv("KIBA_Unique_Proteins.xlsx - Sheet1.csv", engine = 'python')
    # target_proteins = kiba_df["ProteinID"].unique()



    # target_proteins = davis_df["uniprot_id"]
    # target_proteins = target_proteins.unique()
    # target_proteins = list(set(base_proteins) | set(discovered_enzymes))
    # target_proteins = ["P08684"]
    
    # print(f"=== Processing {len(target_proteins)} Primary Proteins (Refreshed List) ===")
    # protein_proc = ProteinProcessor()
    # # Pass 1: Fetch interactions for primary targets
    # protein_proc.process_proteins(target_proteins, fetch_interactions=True)
    
    
    # Save Protein Data
    # for name, df in protein_proc.get_dataframes().items():
    #     print(f"Saving {name}.csv...")
    #     df.to_csv(f"data_kiba/{name}.csv", index=False, mode = 'a')

    # for name, df in protein_proc.get_metadata_dataframes().items():
    #     print(f"Saving {name}.csv...")
    #     df.to_csv(f"data_kiba/{name}.csv", index=False, mode = 'a')

    # Save Drug Data
    # for name, df in drug_proc.get_dataframes().items():
    #     if name == "drug" or name == "drug_side_effects":
    #     # print(f"Saving {name}.csv...")
    #         df.to_csv(f"data_kiba/{name}.csv", index=False, header = False, mode = 'a')

    # print("\nDone! datasets created in 'data_kiba/' directory.")
    '''

# Drug processing below part 3
    # df = pd.read_csv("metz.txt", engine = 'python', sep = '\t')
    # target_drugs = ["CHEMBL941", "CHEMBL1459", "CHEMBL25"]
    # modify the target_drugs
    target_drugs = list(df['chembl_id'].unique())
    total_enzymes = []
    drugProc = DrugProcessor()



    for cid in tqdm(target_drugs, desc = "Processing Drugs: "):
        if not cid:
            continue

        if cid == target_drugs[0]:
            header_condition = True
        else:
            header_condition = False


        data_store = drugProc.process_single_drug(cid)

        if not data_store:
            continue

        for name, df in {name: pd.DataFrame(rows) for name, rows in data_store.items()}.items():
            df.to_csv(f"data/{name}.csv", index=False, header = header_condition, mode = 'a')

        del data_store
        gc.collect()

    discovered_enzymes = drugProc.get_discovered_enzyme_ids()
    for i in discovered_enzymes:    
        total_enzymes.append(i)

    
    print("Enzymes = ", total_enzymes)
    # Calculate Similarity after all drugs processed
    drugProcSimilarity = DrugProcessor()
    data_store = drugProcSimilarity.process_similarity_for_all(target_drugs)
    for name, df in {name: pd.DataFrame(rows) for name, rows in data_store.items()}.items():
        df.to_csv(f"data/{name}.csv", index=False)
    del data_store
    gc.collect()




# Protein processing below part 3
    # target_proteins = ['B7Z8N6', 'D6RGW1', 'Q573B4', 'P36507', 'D7R527', 'D6RG11', 'Q13131', 'P17612', 'P51817', 'O14965', 'M0R3E6', 'A0AA34QVH7', 'A0A0D9SFF9', 'A0A2Q3DQE3', 'A0A0A0MRJ0', 'H3BP07', 'Q16513', 'C9J9P1', 'P05129', 'Q13153', 'H3BLV9', 'H0YNT6', 'P35916', 'D6RAU1', 'Q96KB5', 'O43293', 'Q9Y463', 'Q506Q0', 'A0A0U2VU55', 'A0AAQ5BHR6', 'M0R3G6', 'B1AMW7', 'Q9Y4K4', 'P50750', 'P34947', 'A0A7I2V2J2', 'A0A7I2V3K8', 'A0A0U1RQV4', 'A0A8V8TNA1', 'A0A087WZK4', 'Q8TDC3', 'P42685', 'J3QQI9', 'A8MY48', 'A0A5F9ZH21', 'B1AC84', 'P52333', 'Q16620', 'Q9Y6M4', 'K7EIN7', 'A0A8V8TPW6', 'Q12851', 'Q8N4C8', 'Q6NUK7', 'E9PKQ3', 'A0A2R8YGB7', 'H3BPW1', 'Q9BXA7', 'D3JEN2', 'P11362', 'A0A0D9SFP6', 'C9JYS6', 'A0A7I2V4T2', 'H7BXX9', 'Q5SY34', 'D6RAF9', 'A0A590UJ43', 'A0A6Q8PF65', 'I3L309', 'Q8IZV4', 'Q499G7', 'K7EIM2', 'O94806', 'A0A0S2Z3F9', 'P29376', 'Q02750', 'P11309', 'Q04759', 'A0A7I2V2N0', 'A0A8V8TN50', 'P24941', 'F8VVZ1', 'Q9Y6E0', 'M0QZ82', 'Q9H2G2', 'Q13237', 'P41743', 'A0A2R8YH74', 'F6U4U2', 'I6L9I5', 'F5H6D4', 'A1L4K2', 'A0A6G6D045', 'A0A2R8YCK2', 'Q8NE63', 'Q96GD4', 'A0A7P0Z4D9', 'P08922', 'P49137', 'C3W980', 'A0A090N7W4', 'Q9UEE5', 'Q13689', 'Q9UE13', 'D3DPA4', 'Q9HCP0', 'Q08881', 'P51957', 'Q9H3Y6', 'A0A075B7B4', 'H7C1F0', 'B5BTY5', 'Q07912', 'C9JHS6', 'R4GN93', 'P31749', 'Q86V86', 'B6D4Y6', 'O75116', 'A0A0D9SEY1', 'B7ZL31', 'A0A2R8Y6I6', 'A0ABB0MVC4', 'Q8N6J3', 'E5RII0', 'Q9Y243', 'Q8IU85', 'Q9NYL2', 'Q7KZI7', 'P30530', 'Q8IZL0', 'O00444', 'G3V2D1', 'Q53QS8', 'Q05D26', 'P06241', 'Q71UK5', 'Q13882', 'Q96PF2', 'O15075', 'Q3MS94', 'O96017', 'Q0IJ44', 'Q75MU0', 'Q9H4B4', 'E5RJ77', 'B5BUJ8', 'P10721', 'P08069', 'E7ESA6', 'A0PJF7', 'Q9UEW8', 'G3V105', 'O15111', 'M0R0A7', 'Q9UK79', 'K7EKS5', 'Q96L34', 'Q9H0K1', 'J3KPD6', 'C9JES6', 'A0A0A6YYC0', 'B4DDH2', 'Q9UBE8', 'Q9BWU1', 'O00418', 'E9PQ51', 'Q96RR4', 'A0A286YES9', 'P20815', 'P08684'] #metz
    # target_proteins = ["P00519", "Q9BQB6", "P23219"]
    # modify target_proteins accordingly
    target_proteins = target_proteins + total_enzymes


    proteinProc = ProteinProcessor()


    for uid in tqdm(target_proteins, desc = "Processing Proteins: "):
        if uid == target_proteins[0]:
            header_condition = True
        else:
            header_condition = False


        data_store, metadata_dataframes = proteinProc.process_single_protein(uid, fetch_interactions=True)


        for name, df in {name: pd.DataFrame(rows) for name, rows in data_store.items()}.items():
            df.to_csv(f"data/{name}.csv", index=False, header = header_condition, mode = 'a')

        for name, df in metadata_dataframes.items():
            df.to_csv(f"data/{name}.csv", index=False, header = header_condition,mode = 'a')

        del data_store
        del metadata_dataframes
        gc.collect()

    
# Threading below
    # target_drugs = ["CHEMBL941", "CHEMBL1459", "CHEMBL25"]
    # total_enzymes = []



    # with ThreadPoolExecutor(
    #     max_workers=min(MAX_DRUG_THREADS, len(target_drugs))
    # ) as executor:

    #     futures = {
    #         executor.submit(process_drug_worker, cid, i == 0): cid
    #         for i, cid in enumerate(target_drugs)
    #     }

    #     for future in tqdm(
    #         as_completed(futures),
    #         total=len(futures),
    #         desc="Drugs"
    #     ):
    #         enzymes = future.result()
    #         total_enzymes.extend(enzymes)

    # total_enzymes = list(set(total_enzymes))
    # print("Enzymes: ", total_enzymes)


    # # Similarity
    # drugProcSimilarity = DrugProcessor()
    # data_store = drugProcSimilarity.process_similarity_for_all(target_drugs)
    # for name, df in {name: pd.DataFrame(rows) for name, rows in data_store.items()}.items():
    #     df.to_csv(f"data/{name}.csv", index=False)
    # del data_store
    # gc.collect()



    # # ---------------- PROTEINS ----------------
    # target_proteins = ["P00519", "Q9BQB6", "P23219"]
    # target_proteins += total_enzymes
    # target_proteins = list(dict.fromkeys(target_proteins))



    # with ThreadPoolExecutor(
    #     max_workers=min(MAX_PROTEIN_THREADS, len(target_proteins))
    # ) as executor:

    #     futures = {
    #         executor.submit(process_protein_worker, uid, i == 0): uid
    #         for i, uid in enumerate(target_proteins)
    #     }

    #     for _ in tqdm(
    #         as_completed(futures),
    #         total=len(futures),
    #         desc="Proteins"
    #     ):
    #         pass


    # print("Done")



if __name__ == "__main__":
    main()
