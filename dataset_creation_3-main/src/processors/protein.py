from ..clients import uniprot, pdb, kegg, quickgo, ensembl, stringdb
from tqdm import tqdm
import json
import pandas as pd

class ProteinProcessor:
    def __init__(self):
        self.pathway_lookup = kegg.fetch_kegg_lookup(kegg.KEGG_LIST_PATHWAY, prefix="path:")
        self.disease_lookup = kegg.fetch_kegg_lookup(kegg.KEGG_LIST_DISEASE, prefix="ds:")
        
        self.discovered_interaction_partners = set()

    def process_proteins(self, uniprot_ids, fetch_interactions=True):
        for uid in tqdm(uniprot_ids, desc="Processing Proteins"):
            data = uniprot.fetch_uniprot_data(uid)
            if not data:
                continue

            self._extract_core(data)
            # print("Done 1")
            
            # Features (Motif, Domain, Binding site, Active site)
            # self._extract_features(data, uid)
            # print("Done 2")

            # self._extract_claims(data, uid)
            # print("Done 3")

            # Fetch Variants (Ensembl) & Interactions (STRING) & Tissues (HPA)
            gene_name = data.get("genes", [{}])[0].get("geneName", {}).get("value")
            
            # Get Ensembl ID for HPA
            ensembl_id = None
            for db in data.get("uniProtKBCrossReferences", []):
                if db.get("database") == "Ensembl":
                    # The 'id' field is usually the Transcript ID (ENST).
                    # We need the Gene ID (ENSG) which hides in properties.
                    # properties: [{"key": "protein sequence ID", "value": "..."}, {"key": "gene ID", "value": "ENSG..."}]
                    for prop in db.get("properties", []):
                        if prop.get("key") == "GeneId":
                            val = prop.get("value")
                            if val:
                                ensembl_id = val.split(".")[0] # Strip version
                            break
                    if ensembl_id:
                        break
            
            # self._enrich_variants(uid, gene_name)
            # print("Done 4")

            # self._enrich_variants(uid, gene_name)
            # print("Done 5")
            # if fetch_interactions:
                # self._enrich_string_interactions(uid, gene_name)
                # print("Done 6")
            # self._enrich_tissue_data(uid, ensembl_id)
            # print("Done 7")
            
            pdb_ids, kegg_ids, go_ids = self._extract_external_ids(data)
            
            self._process_pdb(pdb_ids, uid)
            # print("Done 8")
            # self._process_kegg(kegg_ids, uid)
            # print("Done 9")
            # self._process_go(go_ids)
            # print("Done 10")

            # for name, df in self.get_dataframes().items():
            #     if name == "protein_structures":
            #         df.to_csv(f"data_kiba/{name}.csv", index=False)

            # for name, df in self.get_metadata_dataframes().items():
            #     if name == "protein_structures":
            #         df.to_csv(f"data_kiba/{name}.csv", index=False)
            # print("Data written")

    def process_single_protein(self, uid, fetch_interactions=True):
        pathway_metadata_cache = {}
        disease_metadata_cache = {}
        tissue_cache = {} # Name -> UBERON ID
        
        go_term_cache = {}

        data_store = {
            "protein": [],
            "protein_features": [],
            "protein_claims": [],
            "protein_structures": [],
            "protein_kegg_pathways": [],
            "protein_kegg_diseases": [],
            "protein_interactions": [],
            "protein_go": [],
            "ontology_edges": [],
            "protein_tissue": [],
            "protein_variants": [],
            "tissue_nodes": [],
            "tissue_hierarchy": []
        }
        
        

        data = uniprot.fetch_uniprot_data(uid)
        if not data:
            return

        self._extract_core(data, data_store)

        
        # Features (Motif, Domain, Binding site, Active site)
        self._extract_features(data, uid, data_store)


        self._extract_claims(data, uid, data_store)


        # Fetch Variants (Ensembl) & Interactions (STRING) & Tissues (HPA)
        gene_name = data.get("genes", [{}])[0].get("geneName", {}).get("value")
        
        # Get Ensembl ID for HPA
        ensembl_id = None
        for db in data.get("uniProtKBCrossReferences", []):
            if db.get("database") == "Ensembl":
                # The 'id' field is usually the Transcript ID (ENST).
                # We need the Gene ID (ENSG) which hides in properties.
                # properties: [{"key": "protein sequence ID", "value": "..."}, {"key": "gene ID", "value": "ENSG..."}]
                for prop in db.get("properties", []):
                    if prop.get("key") == "GeneId":
                        val = prop.get("value")
                        if val:
                            ensembl_id = val.split(".")[0] # Strip version
                        break
                if ensembl_id:
                    break
        
      
        self._enrich_variants(uid, gene_name, data_store)

        if fetch_interactions:
            self._enrich_string_interactions(uid, gene_name, data_store)

        self._enrich_tissue_data(uid, ensembl_id, data_store, tissue_cache)

        
        pdb_ids, kegg_ids, go_ids = self._extract_external_ids(data, data_store)
        
        self._process_pdb(pdb_ids, uid, data_store)

        self._process_kegg(kegg_ids, uid, data_store, disease_metadata_cache, pathway_metadata_cache)

        # self._process_go(go_ids, data_store, go_term_cache)

    
        return data_store, {
            "ontology_terms": pd.DataFrame(go_term_cache.values()),
            "ontology_pathways": pd.DataFrame(pathway_metadata_cache.values()),
            "ontology_diseases": pd.DataFrame(disease_metadata_cache.values())
        }
    

    def _extract_core(self, data, data_store):
        data_store["protein"].append({
            "protein_id": data.get("primaryAccession"),
            "protein_name": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
            "gene_name": data.get("genes", [{}])[0].get("geneName", {}).get("value"),
            "organism": data.get("organism", {}).get("scientificName"),
            "taxon_id": data.get("organism", {}).get("taxonId"),
            "sequence_length": data.get("sequence", {}).get("length"),
            "annotation_score": data.get("annotationScore"),
            "protein_existence": data.get("proteinExistence"),
            "sequence": data.get("sequence", {}).get("value"),
        })

    def _extract_features(self, data, pid, data_store):
        valid_types = {"Motif", "Domain", "Binding site", "Active site"}
        
        for f in data.get("features", []):
            ftype = f.get("type")
            if ftype in valid_types:
                data_store["protein_features"].append({
                    "protein_id": pid,
                    "type": ftype,
                    "description": f.get("description"),
                    "start": f.get("location", {}).get("start", {}).get("value"),
                    "end": f.get("location", {}).get("end", {}).get("value"),
                })

    def _extract_claims(self, data, pid, data_store):
        for c in data.get("comments", []):
            ctype = c.get("commentType")
            
            # Interactions - DEPRECATED: Replaced by STRING DB
            if ctype == "INTERACTION":
                pass 
                # for inter in c.get("interactions", []):
                #     i2 = inter.get("interactantTwo", {})
                #     self.data_store["protein_interactions"].append({
                #         "protein_id_1": pid,
                #         "protein_id_2": i2.get("uniProtKBAccession"),
                #         "gene_name_2": i2.get("geneName"),
                #         "num_experiments": inter.get("numberOfExperiments")
                #     })
            
            # Generic Text Claims
            for t in c.get("texts", []):
                data_store["protein_claims"].append({
                    "protein_id": pid,
                    "claim_type": ctype,
                    "claim_text": t.get("value"),
                })
                

    def _extract_external_ids(self, data, data_store):
        pdb_ids, kegg_ids, go_ids = set(), set(), set()
        pid = data.get("primaryAccession")
        
        for ref in data.get("uniProtKBCrossReferences", []):
            if ref["database"] == "PDB":
                pdb_ids.add(ref["id"])
            elif ref["database"] == "KEGG":
                kegg_ids.add(ref["id"])
            elif ref["database"] == "GO":
                gid = ref["id"]
                go_ids.add(gid)
                # Store direct link
                data_store["protein_go"].append({
                    "protein_id": pid,
                    "go_id": gid
                })
        return pdb_ids, kegg_ids, go_ids

    def _process_pdb(self, pdb_ids, protein_id, data_store):
        best = None
        for pdb_id in pdb_ids:
            data = pdb.fetch_pdb_entry(pdb_id)
            if not data:
                continue
                
            res = data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
            method = data.get("exptl", [{}])[0].get("method")

            if res is not None:
                if best is None or res < best["resolution"]:
                    # Detailed extraction
                    info = data.get("rcsb_entry_info", {})
                    kw = data.get("struct_keywords", {})
                    refine = data.get("refine", [{}])[0]
                    
                    # Ligands
                    ligands = info.get("nonpolymer_bound_components", [])
                    ligands_str = ", ".join(ligands) if ligands else None
                    
                    best = {
                        "protein_id": protein_id,
                        "pdb_id": pdb_id,
                        "title": data.get("struct", {}).get("title") or data.get("rcsb_primary_citation", {}).get("title"),
                        "method": method,
                        "resolution": res,
                        "structure_source": "PDB",
                        # New Fields
                        "molecular_weight": info.get("molecular_weight"),
                        "complex_type": info.get("selected_polymer_entity_types"), # e.g. Protein/DNA
                        "bound_ligands": ligands_str,
                        "keywords": kw.get("pdbx_keywords"),
                        "r_factor": refine.get("ls_rfactor_rwork")
                    }
        if best:
            data_store["protein_structures"].append(best)

    def _process_kegg(self, kegg_ids, protein_id, data_store, disease_metadata_cache, pathway_metadata_cache):
        for kegg_id in kegg_ids:
            # Pathways
            for pathway_id in kegg.fetch_pathways(kegg_id):
                if pathway_id not in pathway_metadata_cache:
                    name = self.pathway_lookup.get(pathway_id) or self.pathway_lookup.get(pathway_id.replace("path:", ""))
                    # Enrich with description
                    details = kegg.fetch_kegg_entry(pathway_id)
                    description = details.get("DESCRIPTION")
                    
                    pathway_metadata_cache[pathway_id] = {
                        "pathway_id": pathway_id, 
                        "name": name,
                        "description": description
                    }
                
                data_store["protein_kegg_pathways"].append({
                    "protein_id": protein_id,
                    "kegg_gene_id": kegg_id,
                    "pathway_id": pathway_id,
                })

            # Diseases
            for disease_id in kegg.fetch_diseases(kegg_id):
                if disease_id not in disease_metadata_cache:
                    name = self.disease_lookup.get(disease_id) or self.disease_lookup.get(disease_id.replace("ds:", ""))
                    # Enrich with description
                    details = kegg.fetch_kegg_entry(disease_id)
                    description = details.get("DESCRIPTION")

                    disease_metadata_cache[disease_id] = {
                        "disease_id": disease_id, 
                        "name": name,
                        "description": description
                    }

                data_store["protein_kegg_diseases"].append({
                    "protein_id": protein_id,
                    "kegg_gene_id": kegg_id,
                    "disease_id": disease_id,
                })

    def _process_go(self, go_ids, data_store, go_term_cache):
        for go_id in go_ids:
            if go_id not in go_term_cache:
                term = quickgo.fetch_term(go_id)
                if term:
                    go_term_cache[go_id] = {
                        "go_id": go_id,
                        "name": term.get("name"),
                        "aspect": term.get("aspect"),
                        "definition": term.get("definition", {}).get("text"),
                    }
                    parents_data = self._fetch_go_parents_recursive_step(go_id)
                    data_store["ontology_edges"].extend(parents_data)

    def _fetch_go_parents_recursive_step(self, go_id):
        # Implementation of parent fetching that matches main.py logic (single level or using ancestors)
        # main.py actually fetched ancestors, then checked children of ancestors for relation.
        # Here we can simplify or direct import that logic. 
        # For now, I'll basically replicate main.py's specific logic for edge finding.
        
        edges = []
        ancestors = quickgo.fetch_ancestors(go_id)
        if not ancestors:
            return edges
        
        # We need to find which ancestor is a direct parent.
        # This part of main.py was a bit heavy (multi call).
        # We can implement it in quickgo.py or here.
        # Let's rely on quickgo helper.
        
        # Actually logic in main.py:
        # 1. get ancestors
        # 2. get details of ancestors
        # 3. iterate ancestor children to see if CHILD == original GO_ID
        
        ancestors = [a for a in ancestors if a != go_id]
        if not ancestors:
            return edges

        # Fetch details for all ancestors
        terms_data = quickgo.fetch_multiple_terms(ancestors)
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
        return edges

    def get_dataframes(self, data_store):
        import pandas as pd
        result = {
            name: pd.DataFrame(rows) for name, rows in data_store.items()
        }
        return result
    
    def get_metadata_dataframes(self, go_term_cache, pathway_metadata_cache, disease_metadata_cache):
        import pandas as pd
        return {
            "ontology_terms": pd.DataFrame(go_term_cache.values()),
            "ontology_pathways": pd.DataFrame(pathway_metadata_cache.values()),
            "ontology_diseases": pd.DataFrame(disease_metadata_cache.values())
        }

    def _enrich_variants(self, protein_id, gene_name, data_store):
        from ..clients import myvariant  # Lazy import to avoid circular dependency issues if any
        
        if gene_name:
            variants = ensembl.fetch_variants_by_gene(gene_name)
            
            # Step 1: Filter
            filtered_vars = []
            rsids_to_fetch = []
            
            for v in variants:
                desc = v.get("description", "")
                if not desc: continue
                
                desc_lower = desc.lower()
                
                # Filter Logic
                # We relax the description filter to capture ALL potentially pathogenic variants.
                # We strictly rely on Clinical Significance later.
                
                # if "phenotype not specified" in desc_lower: continue  # Removed to catch more
                # if desc_lower.startswith("annotated by"): continue    # Removed to catch more
                
                # ZAP-70 / General Noise Filter:
                # 1. Drop synonymous variants (noise reduction)
                # conseq = v.get("consequence_type", "").lower()
                # if "synonym" in conseq: continue
                # Actually, let's keep synonymous too if they are pathogenic (e.g. splice).
                # The strict "pathogenic" filter will clean them up.
                pass

                # Keep it
                var_obj = {
                    "protein_id": protein_id,
                    "variant_id": v.get("Variation"),
                    "description": desc,
                    "source": v.get("source"),
                    "external_id": v.get("attributes", {}).get("external_id")
                }
                filtered_vars.append(var_obj)
                
                # Collect rsID if available
                rsid = v.get("Variation")
                if rsid and rsid.startswith("rs"):
                    rsids_to_fetch.append(rsid)
            
            # Step 2: Enrich via MyVariant (Batch)
            enrichment_map = {}
            if rsids_to_fetch:
                unique_rsids = list(set(rsids_to_fetch))
                enrichment_map = myvariant.fetch_variants_metadata(unique_rsids)
            
            # Step 3: Merge and Store
            for v_obj in filtered_vars:
                rsid = v_obj["variant_id"]
                clin_sig = None
                
                if rsid in enrichment_map:
                    meta = enrichment_map[rsid]
                    clin_sig = meta.get("clinical_significance")
                    v_obj["clinical_significance"] = clin_sig
                    v_obj["consequence"] = meta.get("consequence")
                    v_obj["protein_change"] = meta.get("protein_change")
                
                # POST-ENRICHMENT FILTER:
                # POST-ENRICHMENT FILTER:
                # STRICT mode: Only keep "Pathogenic" or "Likely pathogenic"
                if not clin_sig:
                    continue
                
                cs_lower = str(clin_sig).lower()
                if "pathogenic" not in cs_lower:
                    continue

                data_store["protein_variants"].append(v_obj)

    def _enrich_string_interactions(self, protein_id, gene_name, data_store):
        # STRING DB Integration
        if not gene_name:
            return

        interactions = stringdb.fetch_interactions([gene_name])
        
        # Collect STRING IDs (Ensembl Protein) to map
        ensembl_ids_to_map = []
        for i in interactions:
            sid = i.get("stringId_B")
            if sid:
                ensembl_ids_to_map.append(sid)
        
        # Batch Map to UniProt
        id_map = uniprot.map_ensembl_to_uniprot(ensembl_ids_to_map)
        
        for i in interactions:
            string_id_b = i.get("stringId_B")
            p2_uniprot = id_map.get(string_id_b) # May be None if mapping failed
            
            if p2_uniprot:
                self.discovered_interaction_partners.add(p2_uniprot)

            data_store["protein_interactions"].append({
                "protein_id_1": protein_id, # Our UniProt ID
                "protein_id_2_string": string_id_b,
                "protein_id_2": p2_uniprot, # NEW: Resolution for ingestor
                "gene_name_2": i.get("preferredName_B"),
                "combined_score": i.get("score"), # 0-1
                "experimental_score": i.get("escore"),
                "database_score": i.get("ascore"),
                "textmining_score": i.get("tscore")
            })

    def get_discovered_interaction_partners(self):
        return list(self.discovered_interaction_partners)

    def _enrich_tissue_data(self, protein_id, ensembl_id, data_store, tissue_cache):
        from ..clients import hpa, ols
        
        # 1. Fetch from HPA
        if not ensembl_id:
             return # No Ensembl ID, can't check HPA reliabley
             
        hpa_data = hpa.fetch_hpa_data(ensembl_id)
        # Result: list of {tissue, value, unit, source}
        
        for item in hpa_data:
            tissue_name = item["tissue"].lower()
            
            # 2. Resolve to Ontology (OLS)
            if tissue_name in tissue_cache:
                term_info = tissue_cache[tissue_name]
            else:
                term_info = ols.search_uberon_by_name(tissue_name)
                if term_info:
                    tissue_cache[tissue_name] = term_info
            
            # If resolved, we proceed. If not, we might skip or store raw name?
            # Let's skip to ensure high quality "Smart" graph.
            if not term_info:
                continue
                
            uberon_id = term_info["id"] # UBERON:1234
            
            # 3. Store Edge (Protein -> Tissue)
            data_store["protein_tissue"].append({
                "protein_id": protein_id,
                "tissue_id": uberon_id,
                "value": item["value"],
                "unit": item["unit"],
                "source": "HPA"
            })
            
            # 4. Store Node (Tissue) - Uniquely?
            # We add to list; ingestor will MERGE. 
            data_store["tissue_nodes"].append({
                "id": uberon_id,
                "name": term_info["name"],
                "description": term_info["description"]
            })
            
            # 5. Build Hierarchy (Tissue -> Parents)
            # Only do this once per tissue (check if we already fetched parents?)
            # Since self.tissue_cache stores just term_info, maybe we can flag if parents fetched?
            # But duplicate fetches in processor are okay, we can rely on caching inside ols if we wanted.
            # For now, let's just fetch.
            
            # Optimization: Check if we have seen this UBERON ID in this run?
            # self.tissue_cache is populated. We can add a 'parents_fetched' flag to it? 
            # Or just fetch. OLS calls are fast enough or volume is low (tissues distinct count < 100).
            if not term_info.get("parents_fetched"):
                 parents = ols.fetch_term_parents(term_info["iri"])
                 for p in parents:
                     # Store Parent Node
                     data_store["tissue_nodes"].append({
                         "id": p["id"],
                         "name": p["name"],
                         "description": p["description"]
                     })
                     # Store Edge (Child -> Parent)
                     data_store["tissue_hierarchy"].append({
                         "source": uberon_id,
                         "target": p["id"],
                         "type": "IS_A"
                     })
                 term_info["parents_fetched"] = True

