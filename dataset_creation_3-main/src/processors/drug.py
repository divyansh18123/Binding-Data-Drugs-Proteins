from ..clients import chembl, openfda, mychem
from tqdm import tqdm
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit import DataStructs

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

class DrugProcessor:
    def __init__(self):
        self.mychem_client = mychem.MyChemClient()
        
        self.discovered_enzyme_ids = set()
        

    def process_drugs(self, chembl_ids):
        for cid in tqdm(chembl_ids, desc="Processing Drugs"):
            mol = chembl.fetch_molecule(cid)
            if not mol:
                continue
            
            self._extract_drug_info(mol, cid)
            # print("Done 1")
            # self._process_targets(cid)
            # print("Done 2")
            # Enrichment with Error Handling
            # If one API fails (e.g. flaky internet), we shouldn't crash the whole pipeline.
            try:
                self._enrich_metabolism(cid)
                # print("Done 3")
            except Exception as e:
                print(f"Warning: Failed to enrich metabolism for {cid}: {e}")

            try:
                self._enrich_side_effects(mol.get('pref_name'), cid)
                # print("Done 4")
            except Exception as e:
                print(f"Warning: Failed to enrich side effects for {cid}: {e}")
            
            try:
                self._enrich_from_mychem(cid, mol.get('pref_name'))
                # print("Done 5")
            except Exception as e:
                print(f"Warning: Failed to enrich from MyChem for {cid}: {e}")
            
            try:
                self._enrich_clinical_labels(mol.get('pref_name'), cid)
                # print("Done 6")
            except Exception as e:
                print(f"Warning: Failed to enrich clinical labels for {cid}: {e}")
            
            # Store mol for similarity
            smiles = mol.get('molecule_structures', {}).get('canonical_smiles')
            if RDKIT_AVAILABLE and smiles:
                m = Chem.MolFromSmiles(smiles)
                if m:
                    self.rdkit_mols[cid] = m
                    
            try:
                self._enrich_atc_hierarchy(cid, mol.get('atc_classifications', []))
                # print("Done 7")
            except Exception as e:
                print(f"Warning: Failed to enrich ATC hierarchy for {cid}: {e}")
            
            try:
                self._enrich_adme(cid)
                # print("Done 8")
            except Exception as e:
                print(f"Warning: Failed to enrich ADME for {cid}: {e}")


            # for name, df in self.get_dataframes().items():
            #     if name == "drug_protein_edges":
            #         df.to_csv(f"data_kiba/{name}.csv", index=False)


        # Calculate Similarity after all drugs processed
        # if RDKIT_AVAILABLE and len(self.rdkit_mols) > 1:
        #     self._calculate_similarity()

        # for name, df in self.get_dataframes().items():
        #     if name == "drug_similarity":
        #         df.to_csv(f"data_kiba/{name}.csv", index=False)
        # print("Data written")


    def process_single_drug(self, cid):
        data_store = {
            "drug": [],
            "drug_protein_edges": [],
            "drug_indication": [],
            "drug_atc": [],
            "drug_metabolism": [],
            "drug_side_effects": [],
            "drug_similarity": [],
            "drug_pharmacogenomics": [],
            "drug_clinical_labels": [],
            "enzyme": [],
            "atc_nodes": [],
            "atc_edges": [],
            "drug_atc_hierarchy": [],
            "adme_nodes": [],
            "drug_adme_edges": []
        }
        target_cache = {}
        rdkit_mols = {} # Store mols for similarity
        


        mol = chembl.fetch_molecule(cid)
        if not mol:
            return
        
        self._extract_drug_info(mol, cid, data_store)

        # self._process_targets(cid, data_store, target_cache)

        # Enrichment with Error Handling
        # If one API fails (e.g. flaky internet), we shouldn't crash the whole pipeline.
        try:
            self._enrich_metabolism(cid, data_store, target_cache)

        except Exception as e:
            print(f"Warning: Failed to enrich metabolism for {cid}: {e}")

        try:
            self._enrich_side_effects(mol.get('pref_name'), cid, data_store)

        except Exception as e:
            print(f"Warning: Failed to enrich side effects for {cid}: {e}")
        
        try:
            self._enrich_from_mychem(cid, data_store, mol.get('pref_name'))

        except Exception as e:
            print(f"Warning: Failed to enrich from MyChem for {cid}: {e}")
        
        try:
            self._enrich_clinical_labels(mol.get('pref_name'), cid, data_store)

        except Exception as e:
            print(f"Warning: Failed to enrich clinical labels for {cid}: {e}")
        
        # Store mol for similarity
        # smiles = mol.get('molecule_structures', {}).get('canonical_smiles')
        # if RDKIT_AVAILABLE and smiles:
        #     m = Chem.MolFromSmiles(smiles)
        #     if m:
        #         self.rdkit_mols[cid] = m
                
        try:
            self._enrich_atc_hierarchy(cid, mol.get('atc_classifications', []), data_store)

        except Exception as e:
            print(f"Warning: Failed to enrich ATC hierarchy for {cid}: {e}")
        
        try:
            self._enrich_adme(cid, data_store)

        except Exception as e:
            print(f"Warning: Failed to enrich ADME for {cid}: {e}")

        return data_store

        
    def process_similarity_for_all(self, target_drugs):
        rdkit_mols = {}
        data_store = {"drug_similarity": []}

        for cid in tqdm(target_drugs, desc = "Processing Similarity: "):
            mol = chembl.fetch_molecule(cid)
            if not mol:
                continue

            smiles = mol.get('molecule_structures', {}).get('canonical_smiles')
            if RDKIT_AVAILABLE and smiles:
                m = Chem.MolFromSmiles(smiles)
                if m:
                    rdkit_mols[cid] = m
                
        if RDKIT_AVAILABLE and len(rdkit_mols) > 1:
            return self._calculate_similarity(data_store, rdkit_mols)
        


    def _extract_drug_info(self, mol, cid, data_store):
        props = mol.get('molecule_properties', {})
        data_store["drug"].append({
            "chembl_id": cid,
            "pref_name": mol.get('pref_name'),
            "molecule_type": mol.get('molecule_type'),
            "mw": props.get('full_mwt'),
            "logp": props.get('alogp'),
            "psa": props.get('psa'),
            "hba": props.get('hba'),
            "hbd": props.get('hbd'),
            "canonical_smiles": mol.get('molecule_structures', {}).get('canonical_smiles'),
        })

        # ATC Codes (Edges to Leaf Nodes) - We now process full hierarchy separately
        # But we still keep the direct link for legacy or simple querying if needed.
        # However, _enrich_atc_hierarchy handles the links. 
        # We can keep 'drug_atc' as a simple link to the L5 code.
        atc_codes = mol.get('atc_classifications', [])
        for code in atc_codes:
            data_store["drug_atc"].append({
                "chembl_id": cid,
                "atc_code": code
            })
            
        # Indications
        indications = chembl.fetch_indications(cid)
        for ind in indications:
            data_store["drug_indication"].append({
                "chembl_id": cid,
                "mesh_id": ind.get('mesh_id'),
                "efo_id": ind.get('efo_id'),
                "indication_term": ind.get('mesh_heading') or ind.get('efo_term'),
                "max_phase": ind.get('max_phase_for_ind')
            })

    def _enrich_atc_hierarchy(self, cid, atc_codes, data_store):
        """
        Fetches full hierarchy for each code and builds the graph segments.
        Nodes: ATC_Level1 -> ... -> ATC_Level4 (Level 5 is the code itself usually)
        Edges: Drug -> Level4 (or 5), Level4 -> Level3 ...
        """
        # We need a set to avoid duplicate node creation if multiple drugs share classes
        # But here we are processing per drug. The ingestor handles MERGE.
        # We just produce the rows.
        
        for code in atc_codes:
            details = chembl.fetch_atc_details(code)
            # details is a list of dicts (usually 1 per code) containing level info
            for d in details:
                # Level 1
                data_store.setdefault("atc_nodes", []).append({
                    "code": d['level1'], "description": d['level1_description'], "level": 1
                })
                # Level 2
                data_store["atc_nodes"].append({
                    "code": d['level2'], "description": d['level2_description'], "level": 2
                })
                # Level 3
                data_store["atc_nodes"].append({
                    "code": d['level3'], "description": d['level3_description'], "level": 3
                })
                # Level 4
                data_store["atc_nodes"].append({
                    "code": d['level4'], "description": d['level4_description'], "level": 4
                })
                
                # Edges: Hierarchy
                data_store.setdefault("atc_edges", []).extend([
                    {"source": d['level2'], "target": d['level1'], "type": "IS_CHILD_OF"},
                    {"source": d['level3'], "target": d['level2'], "type": "IS_CHILD_OF"},
                    {"source": d['level4'], "target": d['level3'], "type": "IS_CHILD_OF"},
                ])
                
                # Edge: Drug -> Level 4 (The most specific class usually available as a group)
                # Note: 'code' (L5) is the specific substance. 
                # If we want to link Drug to its Class, we link to Level 4.
                data_store.setdefault("drug_atc_hierarchy", []).append({
                    "drug_id": cid,
                    "atc_code": d['level4']
                })

    def _enrich_adme(self, cid, data_store):
        adme_data = chembl.fetch_adme_data(cid)
        for act in adme_data:
            param_name = act.get('standard_type') or act.get('type')
            if not param_name: continue
            
            # Create Parameter Node (merged by name)
            data_store.setdefault("adme_nodes", []).append({
                "name": param_name 
            })
            
            # Create Rich Edge
            data_store.setdefault("drug_adme_edges", []).append({
                "drug_id": cid,
                "param_name": param_name,
                "value": act.get('standard_value') or act.get('value'),
                "unit": act.get('standard_units') or act.get('units'),
                "organism": act.get('target_organism') or "Unknown",
                "assay_type": "ADME"
            })

    def _process_targets(self, cid, data_store, target_cache):
        # 1. Try Mechanisms first
        mechanisms = chembl.fetch_mechanisms(cid)
        found_targets = False
        
        for mech in mechanisms:
            target_chembl_id = mech.get('target_chembl_id')
            if target_chembl_id:
                self._add_edge(cid, target_chembl_id, "mechanism", mech.get('mechanism_of_action'), data_store, target_cache)
                found_targets = True
        
        # 2. If no mechanisms (or valid targets from them), try activities
        if not found_targets: # Or maybe always? User asked to "figure out... relations". Let's do both but prioritize mechanism.
            # Let's add activities as well, but maybe filter for strong binding.
            activities = chembl.fetch_activities(cid)
            # Filter for human targets and reasonable types
            valid_types = ['IC50', 'Ki', 'EC50', 'Kd']
            
            for act in activities:
                # Basic filter: must have a target and a standard type we care about
                if act.get('standard_type') in valid_types and act.get('target_chembl_id'):
                    # We could also check pChEMBL value if we want strict quality
                    self._add_edge(
                        cid, 
                        act.get('target_chembl_id'), 
                        "activity", 
                        f"{act.get('standard_type')} = {act.get('standard_value')} {act.get('standard_units')}",
                        data_store,
                        target_cache
                    )

    def _add_edge(self, drug_id, target_chembl_id, relation_type, relation_desc, data_store, target_cache):
        # Convert Target ChEMBL ID to Uniprot
        uniprot_ids = self._get_uniprots_for_target(target_chembl_id, target_cache)
        
        for uid in uniprot_ids:
            edge_data = {
                "drug_chembl_id": drug_id,
                "protein_uniprot_id": uid,
                "target_chembl_id": target_chembl_id,
                "relation_type": relation_type, # mechanism or activity
                "description": relation_desc,
            }
            # Add structured affinity data if available
            data_store["drug_protein_edges"].append(edge_data)

    def _get_uniprots_for_target(self, target_chembl_id, target_cache):
        if target_chembl_id in target_cache:
            return target_cache[target_chembl_id]
        
        target_info = chembl.fetch_target(target_chembl_id)
        uniprots = set()
        
        if target_info:
            # Check organism
            if target_info.get('organism') == 'Homo sapiens':
                for comp in target_info.get('target_components', []):
                    acc = comp.get('accession')
                    if acc:
                        uniprots.add(acc)
        
        target_cache[target_chembl_id] = list(uniprots)
        return list(uniprots)

    def _enrich_metabolism(self, cid, data_store, target_cache):
        met_data = chembl.fetch_metabolism(cid)
        for m in met_data:
            enzyme_chembl_id = m.get("target_chembl_id")
            uniprot_id = None
            
            # Logic for Enzyme Nodes
            if enzyme_chembl_id:
                # 1. Resolve to UniProt
                uniprot_ids = self._get_uniprots_for_target(enzyme_chembl_id, target_cache)
                # We usually take the first one if multiple, or create multiple enzyme entries?
                # ChEMBL target = Single Protein usually maps to 1 UniProt.
                uniprot_id = uniprot_ids[0] if uniprot_ids else None
                
                if uniprot_id:
                    self.discovered_enzyme_ids.add(uniprot_id)
            
            data_store["drug_metabolism"].append({
                "drug_chembl_id": cid,
                "enzyme_name": m.get("enzyme_name"),
                "enzyme_chembl_id": enzyme_chembl_id,
                "enzyme_uniprot_id": uniprot_id, # Can be None if lookup failed
                "metabolite_name": m.get("metabolite_name"),
                "mapping_type": "Metabolism" 
            })

    def get_discovered_enzyme_ids(self):
        return list(self.discovered_enzyme_ids)

    def _enrich_from_mychem(self, cid, data_store, drug_name=None):
        """
        Fetches data from MyChem.info and updates:
        1. Drug node properties (PubChem ID, DrugBank ID) - by modifying the last added drug entry.
        2. Pharmacogenomic edges (Drug-Gene).
        """
        data = self.mychem_client.fetch_data(cid)
        
        # Fallback to Name Search if ID failed
        if not data and drug_name:
             print(f"  MyChem ID lookup failed for {cid}, trying name: {drug_name}...")
             data = self.mychem_client.fetch_data_by_name(drug_name)

        if not data:
            return



        # --- Information Extraction ---
        # 1. External IDs
        # JSON Path: hit['chembl']['molecule_chembl_id'] (verified match)
        # JSON Path: hit['pubchem']['cid']
        # JSON Path: hit['drugbank']['id']
        try:
            pubchem_id = data.get('pubchem', {}).get('cid')
        except:
            try:
                if isinstance(data.get('pubchem', {}), list): 
                    pubchem_id = data.get('pubchem', {})[0].get('cid') # Handle list case

            except:
                raise Exception("Could not enrich mychem")
        
        try:
            drugbank_id = data.get('drugbank', {}).get('id')
        except:
            try:
                if isinstance(data.get('drugbank', {}), list): 
                    drugbank_id = data.get('drugbank', {})[0].get('id')

            except:
                raise Exception("Could not enrich mychem")

        # Update the last added drug entry with these new IDs
        # We assume the loop order is maintained and we just added this 'cid'
        if data_store["drug"] and data_store["drug"][-1]["chembl_id"] == cid:
            data_store["drug"][-1]["pubchem_id"] = pubchem_id
            data_store["drug"][-1]["drugbank_id"] = drugbank_id
            
            # Additional enrichment: SIDER side effect count
            # JSON Path: hit['sider'] (List of dicts)
            sider_data = data.get('sider', [])
            if sider_data:
                
                # Extract SIDER side effects
                for s in sider_data:
                    side_effect_name = s.get('side_effect', {}).get('name')
                    if not side_effect_name: continue

                    data_store["drug"][-1]["side_effect_count"] += 1
                    
                    # Extract Rich Attributes
                    frequency = s.get('side_effect', {}).get('frequency')
                    
                    # Placebo is often a boolean in SIDER/MyChem (true = observed in placebo?)
                    # Or sometimes a frequency string. Let's capture what is there.
                    placebo_val = s.get('placebo')
                    if placebo_val is None:
                         # Check nested side_effect.placebo if it exists
                         placebo_val = s.get('side_effect', {}).get('placebo')
                    
                    # MedDRA ID hierarchy: preferred > concept > umls
                    meddra_obj = s.get('meddra', {})
                    meddra_code = meddra_obj.get('preferred_term_code') or meddra_obj.get('concept_id') or meddra_obj.get('umls_id')
                    meddra_type = meddra_obj.get('type')
                    
                    # Normalize Frequency (it can be a string or localized dict)
                    freq_val = frequency
                    if isinstance(frequency, dict):
                        # rare case in SIDER json sometimes
                        freq_val = frequency.get('value') or str(frequency)
                        
                    data_store["drug_side_effects"].append({
                        "drug_chembl_id": cid,
                        "side_effect": side_effect_name,
                        "frequency": freq_val,
                        "placebo_frequency": str(placebo_val) if placebo_val is not None else None,
                        "meddra_code": meddra_code,
                        "meddra_type": meddra_type,
                        "source": "SIDER"
                    })

        # 2. Pharmacogenomics (PharmGKB)
        # JSON Path: hit['pharmgkb'] (Dict or List)
        pg_data = data.get('pharmgkb')
        if pg_data:
            # PharmGKB field can be a single dict or list of dicts. Normalize to list.
            if isinstance(pg_data, dict):
                 pg_data = [pg_data]
            
            for item in pg_data:
                # JSON Path: item['gene'] (List of gene symbols, e.g. ['CYP2D6'])
                genes = item.get('gene', [])
                if isinstance(genes, str): genes = [genes]
                
                # JSON Path: item['name'] (Name of the interaction/variant e.g. "rs1065852")
                # JSON Path: item['association'] (The effect, e.g. "associated with increased response")
                phenotype = item.get('association') or item.get('name')
                
                for gene_symbol in genes:
                    data_store["drug_pharmacogenomics"].append({
                        "drug_id": cid,
                        "gene_symbol": gene_symbol,
                        "phenotype": phenotype,
                        "source": "PharmGKB"
                    })

    def _enrich_side_effects(self, drug_name, cid, data_store):
        effects = openfda.fetch_side_effects(drug_name)
        data_store["drug"][-1]["side_effect_count"] = 0
        for effect in effects:
            data_store["drug"][-1]["side_effect_count"] += 1
            data_store["drug_side_effects"].append({
                "drug_chembl_id": cid,
                "side_effect": effect,
                "source": "OpenFDA"
            })

    def _enrich_clinical_labels(self, drug_name, cid, data_store):
        """
        Fetches full text clinical info (Mechanism, Interactions) from OpenFDA.
        """
        label_data = openfda.fetch_drug_label(drug_name)
        if label_data:
            data_store["drug_clinical_labels"].append({
                "chembl_id": cid,
                "mechanism_of_action": label_data.get("mechanism_of_action"),
                "drug_interactions": label_data.get("drug_interactions"),
                "contraindications": label_data.get("contraindications"),
                "boxed_warning": label_data.get("boxed_warning"),
                "clinical_pharmacology": label_data.get("clinical_pharmacology")
            })



    def _calculate_similarity(self, data_store, rdkit_mols):
        # # Fingerprints
        # fps = {cid: AllChem.GetMorganFingerprintAsBitVect(m, 2) for cid, m in self.rdkit_mols.items()}

        # Create Morgan fingerprint generator (radius=2, default fpSize=2048)
        morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)

        # Generate fingerprints
        fps = {
            cid: morgan_gen.GetFingerprint(mol)
            for cid, mol in rdkit_mols.items()
        }

        cids = list(fps.keys())
        
        for i in range(len(cids)):
            for j in range(i + 1, len(cids)):
                id1 = cids[i]
                id2 = cids[j]
                sim = DataStructs.TanimotoSimilarity(fps[id1], fps[id2])
                
                # Store extensive matrix or just threshold? 
                # Storing all pairs for small dataset is fine.
                if sim > 0.0:  # Threshold > 0 to save space if needed
                    data_store["drug_similarity"].append({
                        "drug_id_1": id1,
                        "drug_id_2": id2,
                        "similarity_score": sim
                    })
        return data_store

    def get_dataframes(self, data_store):
        import pandas as pd
        return {
            name: pd.DataFrame(rows) for name, rows in data_store.items()
        }
