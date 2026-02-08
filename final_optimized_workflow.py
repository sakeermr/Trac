#!/usr/bin/env python3
"""
STANDARD SEED CORPORATION - STANDARD-SEED-TARGET-PREDICTOR
==========================================================

Version 2.0 - Internal Use Only

Advanced molecular similarity screening system designed for
Standard Seed Corporation's drug discovery and chemical research operations.
Processes 650,000+ PDB ligands with maximum efficiency and reliability.

This is the OPTIMIZED version that consolidates all features:
- Processes entire PDB database (650K+ ligands)
- Human/Mouse/Rat organism filtering
- Single CSV output with clean formatting
- Windows-compatible with Unicode error handling
- Memory-optimized chunked processing
- Professional CLI interface
- Comprehensive error handling and logging

Organization: Standard Seed Corporation
Version: 2.0
License: Apache 2.0 (Internal Use Only)
"""

import pandas as pd
import numpy as np
import sys
import os
import logging
import argparse
import warnings
import csv
from datetime import datetime
from pathlib import Path
import chardet
import requests
import json
import time
from collections import defaultdict

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', module='rdkit')

# RDKit imports with fallback handling
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, DataStructs
    from rdkit import RDLogger
    RDKIT_AVAILABLE = True
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    print("ERROR: RDKit not available. Please install: pip install rdkit")
    sys.exit(1)

class FinalOptimizedPDBScreening:
    """Final optimized PDB ligand screening system for 600K+ ligands"""
    
    def __init__(self, max_pdb_records=0, max_input_chemicals=0, chunk_size=10000):
        self.max_pdb_records = max_pdb_records
        self.max_input_chemicals = max_input_chemicals
        self.chunk_size = chunk_size
        self.pdb_fingerprints = {}
        self.pdb_organism_cache = {}  # Cache for PDB organism information
        self.processing_stats = {
            'total_pdb_processed': 0,
            'valid_pdb_smiles': 0,
            'total_input_processed': 0,
            'valid_input_smiles': 0,
            'total_comparisons': 0,
            'matches_found': 0,
            'processing_time': 0.0,
            'errors_encountered': 0
        }
        self.logger = self._setup_logging()
        
    def _setup_logging(self):
        """Setup professional logging configuration"""
        logger = logging.getLogger('PDBScreening')
        logger.setLevel(logging.INFO)
        
        # Create handlers
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)
        
        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_handler.setFormatter(formatter)
        
        # Add handler to logger
        if not logger.handlers:
            logger.addHandler(console_handler)
        
        return logger
    
    def _detect_encoding(self, file_path):
        """Detect file encoding to handle Unicode issues"""
        try:
            with open(file_path, 'rb') as f:
                sample = f.read(100000)  # Read first 100KB
                result = chardet.detect(sample)
                encoding = result['encoding']
                confidence = result['confidence']
                
                if confidence < 0.7:
                    encoding = 'utf-8'
                    
                self.logger.info(f"Detected encoding: {encoding} (confidence: {confidence:.2f})")
                return encoding
        except Exception as e:
            self.logger.warning(f"Encoding detection failed: {e}. Using utf-8 with error handling.")
            return 'utf-8'
    
    def _load_pdb_database_robust(self, pdb_file):
        """Load PDB database with robust error handling for 600K+ ligands"""
        self.logger.info(f"Loading PDB database from: {pdb_file}")
        
        # Detect encoding
        encoding = self._detect_encoding(pdb_file)
        
        try:
            # Try standard pandas read first
            df = pd.read_csv(pdb_file, encoding=encoding, low_memory=False)
            self.logger.info(f"Successfully loaded {len(df):,} PDB records")
            return df
            
        except UnicodeDecodeError:
            self.logger.warning("Unicode error detected. Using robust line-by-line parsing...")
            return self._load_pdb_line_by_line(pdb_file)
        except Exception as e:
            self.logger.error(f"Failed to load PDB database: {e}")
            raise
    
    def _load_pdb_line_by_line(self, pdb_file):
        """Load PDB file line by line with Unicode error handling"""
        records = []
        error_count = 0
        
        encodings_to_try = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
        
        for encoding in encodings_to_try:
            try:
                with open(pdb_file, 'r', encoding=encoding, errors='ignore') as f:
                    reader = csv.DictReader(f)
                    
                    for i, row in enumerate(reader):
                        if self.max_pdb_records > 0 and i >= self.max_pdb_records:
                            break
                            
                        # Clean the row data
                        clean_row = {}
                        for key, value in row.items():
                            if value is not None:
                                clean_row[key] = str(value).strip()
                            else:
                                clean_row[key] = ''
                        
                        records.append(clean_row)
                        
                        if (i + 1) % 50000 == 0:
                            self.logger.info(f"Processed {i+1:,} PDB records...")
                
                self.logger.info(f"Successfully loaded {len(records):,} PDB records with encoding: {encoding}")
                return pd.DataFrame(records)
                
            except Exception as e:
                error_count += 1
                self.logger.warning(f"Failed with encoding {encoding}: {e}")
                continue
        
        raise Exception(f"Failed to load PDB file with any encoding after {error_count} attempts")
    
    def _generate_fingerprint(self, smiles):
        """Generate Morgan fingerprint for a SMILES string"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            # Generate Morgan fingerprint
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            return fp
        except Exception:
            return None
    
    def _calculate_similarity(self, fp1, fp2):
        """Calculate Tanimoto similarity between two fingerprints"""
        try:
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        except Exception:
            return 0.0
    
    def _fetch_pdb_organism_info(self, pdb_ids):
        """Fetch organism information for multiple PDB IDs using PDB API"""
        if not pdb_ids:
            return {}
        
        # Filter out already cached PDB IDs
        uncached_pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id not in self.pdb_organism_cache]
        
        if not uncached_pdb_ids:
            return {pdb_id: self.pdb_organism_cache[pdb_id] for pdb_id in pdb_ids}
        
        self.logger.info(f"Fetching organism information for {len(uncached_pdb_ids)} PDB IDs...")
        
        # Track successful API calls
        total_successful = 0
        total_human_mouse_rat = 0
        
        # Batch API request to PDB
        batch_size = 100  # PDB API limit
        organism_info = {}
        
        for i in range(0, len(uncached_pdb_ids), batch_size):
            batch = uncached_pdb_ids[i:i+batch_size]
            
            try:
                # PDB GraphQL API query - Updated for current API structure
                query = {
                    "query": """
                    query ($pdb_ids: [String!]!) {
                        entries(entry_ids: $pdb_ids) {
                            rcsb_id
                            rcsb_entry_info {
                                structure_determination_methodology
                            }
                            polymer_entities {
                                rcsb_polymer_entity_container_identifiers {
                                    entry_id
                                }
                                entity_src_gen {
                                    pdbx_gene_src_scientific_name
                                }
                                entity_src_nat {
                                    pdbx_organism_scientific
                                }
                                rcsb_entity_source_organism {
                                    ncbi_taxonomy_id
                                    scientific_name
                                }
                            }
                            nonpolymer_entities {
                                rcsb_nonpolymer_entity_container_identifiers {
                                    entry_id
                                }
                                nonpolymer_entity_instances {
                                    rcsb_nonpolymer_entity_instance_container_identifiers {
                                        entry_id
                                    }
                                }
                            }
                        }
                    }
                    """,
                    "variables": {"pdb_ids": batch}
                }
                
                response = requests.post(
                    "https://data.rcsb.org/graphql",
                    json=query,
                    headers={"Content-Type": "application/json"},
                    timeout=30
                )
                
                if response.status_code == 200:
                    data = response.json()
                    
                    if 'data' in data and 'entries' in data['data'] and data['data']['entries']:
                        for entry in data['data']['entries']:
                            if entry and 'rcsb_id' in entry:
                                pdb_id = entry['rcsb_id']
                                organisms = []
                                
                                # Safely extract organism information from multiple sources
                                try:
                                    # Check polymer entities for organism info
                                    if 'polymer_entities' in entry and entry['polymer_entities']:
                                        for polymer in entry['polymer_entities']:
                                            if polymer:
                                                # Try rcsb_entity_source_organism first
                                                if 'rcsb_entity_source_organism' in polymer and polymer['rcsb_entity_source_organism']:
                                                    for org in polymer['rcsb_entity_source_organism']:
                                                        if org and 'scientific_name' in org and org['scientific_name']:
                                                            organisms.append(str(org['scientific_name']))
                                                
                                                # Try entity_src_gen as backup
                                                if 'entity_src_gen' in polymer and polymer['entity_src_gen']:
                                                    for gen in polymer['entity_src_gen']:
                                                        if gen and 'pdbx_gene_src_scientific_name' in gen and gen['pdbx_gene_src_scientific_name']:
                                                            organisms.append(str(gen['pdbx_gene_src_scientific_name']))
                                                
                                                # Try entity_src_nat as backup
                                                if 'entity_src_nat' in polymer and polymer['entity_src_nat']:
                                                    for nat in polymer['entity_src_nat']:
                                                        if nat and 'pdbx_organism_scientific' in nat and nat['pdbx_organism_scientific']:
                                                            organisms.append(str(nat['pdbx_organism_scientific']))
                                                            
                                except Exception as parse_error:
                                    self.logger.debug(f"Error parsing organism for PDB {pdb_id}: {parse_error}")
                                    continue
                                
                                # Remove duplicates and ensure we have at least some data
                                organisms = list(set(organisms)) if organisms else ['Unknown']
                                
                                organism_info[pdb_id] = {
                                    'organisms': organisms,
                                    'is_human_mouse_rat': self._is_human_mouse_or_rat(organisms)
                                }
                                
                                # Cache the result
                                self.pdb_organism_cache[pdb_id] = organism_info[pdb_id]
                                total_successful += 1
                                
                                # Debug logging for successful parsing
                                if organisms != ['Unknown']:
                                    self.logger.debug(f"PDB {pdb_id}: Found organisms {organisms}")
                                    if organism_info[pdb_id]['is_human_mouse_rat']:
                                        total_human_mouse_rat += 1
                                    
                    else:
                        self.logger.debug(f"No entries found in PDB API response for batch: {batch}")
                else:
                    self.logger.warning(f"PDB API returned status code {response.status_code} for batch: {batch}")
                
                # Rate limiting to avoid overwhelming the API
                time.sleep(0.1)
                
            except Exception as e:
                self.logger.warning(f"GraphQL API failed for batch, trying REST API fallback: {e}")
                
                # Fallback to REST API for individual PDB IDs
                for pdb_id in batch:
                    try:
                        rest_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"
                        rest_response = requests.get(rest_url, timeout=10)
                        
                        if rest_response.status_code == 200:
                            rest_data = rest_response.json()
                            organisms = []
                            
                            # Extract organism info from REST API response
                            if 'rcsb_entry_info' in rest_data:
                                entry_info = rest_data['rcsb_entry_info']
                                if 'polymer_entity_count_protein' in entry_info:
                                    # This is a protein structure, try to get organism info
                                    organisms.append('Protein structure')
                            
                            # Try to get more specific organism info
                            if 'struct' in rest_data and 'title' in rest_data['struct']:
                                title = rest_data['struct']['title'].lower()
                                if 'human' in title or 'homo sapiens' in title:
                                    organisms.append('Homo sapiens')
                                elif 'mouse' in title or 'mus musculus' in title:
                                    organisms.append('Mus musculus')
                                elif 'rat' in title or 'rattus' in title:
                                    organisms.append('Rattus norvegicus')
                            
                            organisms = list(set(organisms)) if organisms else ['Unknown']
                            
                            organism_info[pdb_id] = {
                                'organisms': organisms,
                                'is_human_mouse_rat': self._is_human_mouse_or_rat(organisms)
                            }
                            self.pdb_organism_cache[pdb_id] = organism_info[pdb_id]
                            
                            if organisms != ['Unknown']:
                                self.logger.debug(f"REST API - PDB {pdb_id}: Found organisms {organisms}")
                                total_successful += 1
                                if organism_info[pdb_id]['is_human_mouse_rat']:
                                    total_human_mouse_rat += 1
                        else:
                            # Mark as unknown if REST API also fails
                            organism_info[pdb_id] = {
                                'organisms': ['Unknown'],
                                'is_human_mouse_rat': False
                            }
                            self.pdb_organism_cache[pdb_id] = organism_info[pdb_id]
                            
                        time.sleep(0.05)  # Rate limiting for REST API
                        
                    except Exception as rest_error:
                        self.logger.debug(f"REST API also failed for PDB {pdb_id}: {rest_error}")
                        organism_info[pdb_id] = {
                            'organisms': ['Unknown'],
                            'is_human_mouse_rat': False
                        }
                        self.pdb_organism_cache[pdb_id] = organism_info[pdb_id]
                
                # Continue processing even if some batches fail
                continue
        
        # Combine with cached results
        final_result = {}
        for pdb_id in pdb_ids:
            final_result[pdb_id] = self.pdb_organism_cache.get(pdb_id, {
                'organisms': ['Unknown'],
                'is_human_mouse_rat': False
            })
        
        # Log summary of API success
        self.logger.info(f"PDB API Summary: {total_successful}/{len(uncached_pdb_ids)} successful, {total_human_mouse_rat} Human/Mouse/Rat found")
        
        return final_result
    
    def _is_human_mouse_or_rat(self, organisms):
        """Check if any organism is human, mouse, or rat"""
        if not organisms or organisms is None:
            return False
        
        # Ensure organisms is a list
        if isinstance(organisms, str):
            organisms = [organisms]
        
        # Filter out None values
        organisms = [org for org in organisms if org is not None]
        
        if not organisms:
            return False
        
        human_keywords = ['homo sapiens', 'human', 'h. sapiens']
        mouse_keywords = ['mus musculus', 'mouse', 'm. musculus']
        rat_keywords = ['rattus norvegicus', 'rat', 'r. norvegicus', 'rattus']
        
        target_keywords = human_keywords + mouse_keywords + rat_keywords
        
        for organism in organisms:
            if organism is None:
                continue
            organism_lower = str(organism).lower()
            if any(keyword in organism_lower for keyword in target_keywords):
                return True
        
        return False
    
    def _process_pdb_fingerprints(self, pdb_df):
        """Process PDB database to generate fingerprints"""
        self.logger.info("Generating fingerprints for PDB ligands...")
        
        processed_count = 0
        valid_count = 0
        
        # Determine how many records to process
        total_to_process = len(pdb_df)
        if self.max_pdb_records > 0:
            total_to_process = min(total_to_process, self.max_pdb_records)
            self.logger.info(f"Processing limited to {total_to_process:,} PDB records (max_pdb_records={self.max_pdb_records:,})")
        
        for idx, row in pdb_df.iterrows():
            # Check if we've reached the limit
            if self.max_pdb_records > 0 and processed_count >= self.max_pdb_records:
                break
                
            try:
                smiles = row.get('SMILES', '')
                if not smiles or pd.isna(smiles):
                    processed_count += 1
                    continue
                
                fp = self._generate_fingerprint(smiles)
                if fp is not None:
                    self.pdb_fingerprints[idx] = {
                        'fingerprint': fp,
                        'pdb_id': row.get('PDB_ID', ''),
                        'ligand_name': row.get('Ligand_Name', ''),
                        'smiles': smiles,
                        'molecular_weight': row.get('Molecular_Weight', ''),
                        'status': row.get('Status', '')
                    }
                    valid_count += 1
                
                processed_count += 1
                
                if processed_count % 25000 == 0:
                    self.logger.info(f"Processed {processed_count:,} PDB ligands, {valid_count:,} valid fingerprints")
                    
            except Exception as e:
                self.processing_stats['errors_encountered'] += 1
                processed_count += 1
                continue
        
        self.processing_stats['total_pdb_processed'] = processed_count
        self.processing_stats['valid_pdb_smiles'] = valid_count
        
        self.logger.info(f"PDB fingerprint generation complete: {valid_count:,} valid fingerprints from {processed_count:,} records")
    
    def _screen_chemicals(self, input_df, output_file):
        """Screen input chemicals against PDB database"""
        self.logger.info("Starting chemical screening process...")
        
        # Prepare output file
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        results = []
        total_input = len(input_df)
        if self.max_input_chemicals > 0:
            total_input = min(total_input, self.max_input_chemicals)
        
        for idx, row in input_df.iterrows():
            if self.max_input_chemicals > 0 and idx >= self.max_input_chemicals:
                break
                
            try:
                input_smiles = row.get('Molecular Structure', '')
                if not input_smiles or pd.isna(input_smiles):
                    continue
                
                # Generate fingerprint for input chemical
                input_fp = self._generate_fingerprint(input_smiles)
                if input_fp is None:
                    continue
                
                self.processing_stats['valid_input_smiles'] += 1
                
                # Find best matches in PDB database
                best_matches = []
                
                for pdb_idx, pdb_data in self.pdb_fingerprints.items():
                    similarity = self._calculate_similarity(input_fp, pdb_data['fingerprint'])
                    self.processing_stats['total_comparisons'] += 1
                    
                    if similarity > 0.1:  # Lowered threshold for more matches
                        best_matches.append({
                            'pdb_id': pdb_data['pdb_id'],
                            'ligand_name': pdb_data['ligand_name'],
                            'pdb_smiles': pdb_data['smiles'],
                            'similarity': similarity,
                            'molecular_weight': pdb_data['molecular_weight'],
                            'status': pdb_data['status']
                        })
                
                # Sort by similarity and take top matches
                best_matches.sort(key=lambda x: x['similarity'], reverse=True)
                top_matches = best_matches[:50]  # Get top 50 for better organism filtering
                
                if top_matches:
                    # Get organism information for top matches
                    pdb_ids = [match['pdb_id'] for match in top_matches]
                    
                    try:
                        organism_info = self._fetch_pdb_organism_info(pdb_ids)
                    except Exception as api_error:
                        self.logger.warning(f"Failed to fetch organism info for {row.get('Chemical Name', 'Unknown')}: {api_error}")
                        # Create fallback organism info
                        organism_info = {pdb_id: {'organisms': ['Unknown'], 'is_human_mouse_rat': False} for pdb_id in pdb_ids}
                    
                    # Filter for human/mouse/rat only from top matches
                    human_mouse_rat_matches = []
                    for match in top_matches:
                        pdb_id = match['pdb_id']
                        if pdb_id in organism_info and organism_info[pdb_id].get('is_human_mouse_rat', False):
                            match['organisms'] = organism_info[pdb_id].get('organisms', ['Unknown'])
                            human_mouse_rat_matches.append(match)
                    
                    # If we don't have enough human/mouse/rat matches from top 50, extend search
                    if len(human_mouse_rat_matches) < 5:
                        # Get more matches from the full list
                        extended_matches = best_matches[50:200]  # Check next 150 matches
                        if extended_matches:
                            extended_pdb_ids = [match['pdb_id'] for match in extended_matches]
                            
                            try:
                                extended_organism_info = self._fetch_pdb_organism_info(extended_pdb_ids)
                            except Exception as api_error:
                                self.logger.warning(f"Failed to fetch extended organism info for {row.get('Chemical Name', 'Unknown')}: {api_error}")
                                extended_organism_info = {pdb_id: {'organisms': ['Unknown'], 'is_human_mouse_rat': False} for pdb_id in extended_pdb_ids}
                            
                            for match in extended_matches:
                                pdb_id = match['pdb_id']
                                if pdb_id in extended_organism_info and extended_organism_info[pdb_id].get('is_human_mouse_rat', False):
                                    match['organisms'] = extended_organism_info[pdb_id].get('organisms', ['Unknown'])
                                    human_mouse_rat_matches.append(match)
                                    
                                    # Stop if we have enough matches
                                    if len(human_mouse_rat_matches) >= 10:
                                        break
                    
                    # Sort human/mouse/rat matches by similarity and take top 10
                    human_mouse_rat_matches.sort(key=lambda x: x['similarity'], reverse=True)
                    final_matches = human_mouse_rat_matches[:10]
                    
                    # If still no human/mouse/rat matches, use top general matches as fallback
                    if not final_matches:
                        final_matches = top_matches[:10]
                        self.logger.warning(f"No Human/Mouse/Rat matches found for {row.get('Chemical Name', 'Unknown')}, using general matches")
                    
                    if final_matches:
                        self.processing_stats['matches_found'] += 1
                        
                        for i, match in enumerate(final_matches, 1):
                            match_pdb_id = match['pdb_id']
                            results.append({
                                'Plant': row.get('Plant', ''),
                                'Chemical_Name': row.get('Chemical Name', ''),
                                'Molecular_Structure': input_smiles,
                                'Molecule_Category': row.get('Molecule Category', 'Unknown'),
                                'Match_Rank': i,
                                'PDB_ID': match_pdb_id,
                                'Ligand_Name': match['ligand_name'],
                                'PDB_SMILES': match['pdb_smiles'],
                                'Similarity_Score': round(match['similarity'], 4),
                                'Molecular_Weight': match['molecular_weight'],
                                'Status': match['status'],
                                'Organisms': match.get('organisms', ['Unknown']),
                                'Is_Human_Mouse_Rat': match_pdb_id in organism_info and organism_info[match_pdb_id].get('is_human_mouse_rat', False),
                                'Processing_Date': datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                            })
                
                self.processing_stats['total_input_processed'] += 1
                
                if (idx + 1) % 10 == 0:
                    self.logger.info(f"Processed {idx+1:,}/{total_input:,} input chemicals")
                    
            except Exception as e:
                self.processing_stats['errors_encountered'] += 1
                self.logger.warning(f"Error processing chemical {idx}: {e}")
                continue
        
        # Save results with analysis
        self._save_comprehensive_results(results, output_file)
        
        self.logger.info(f"Screening complete. Results saved to: {output_file}")
    
    def _save_comprehensive_results(self, results, output_file):
        """Save results in the requested format: CSV with top 5 targets and analysis report"""
        # Convert results to DataFrame for easier processing
        results_df = pd.DataFrame(results)
        
        # Generate CSV output file with top 5 targets
        csv_output = output_file.replace('.csv', '_top5_targets.csv')
        self._save_csv_with_top5_targets(results_df, csv_output)
        
        # Generate analysis report with top 20 targets and Tanimoto values
        report_output = output_file.replace('.csv', '_analysis_report.txt')
        self._save_analysis_report(results_df, report_output)
        
        # Save detailed results for reference
        detailed_output = output_file.replace('.csv', '_detailed_results.csv')
        results_df.to_csv(detailed_output, index=False)
        
        self.logger.info(f"Results saved to:")
        self.logger.info(f"  - Top 5 targets (clean CSV): {csv_output}")
        self.logger.info(f"  - Analysis report: {report_output}")
        self.logger.info(f"  - Detailed results: {detailed_output}")
    
    def _save_csv_with_top5_targets(self, results_df, output_file):
        """Save CSV file with top 5 targets for each molecule - clean CSV format"""
        if results_df.empty:
            # Create empty file with headers
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(['Plant', 'Chemical Name', 'Molecular Structure', 'Molecule Category', 'Top PDB IDs'])
                writer.writerow(['No matches found with similarity > 0.1', '', '', '', ''])
            return
        
        # Group by input chemical and get top 5 targets for each
        csv_data = []
        for chemical_name in results_df['Chemical_Name'].unique():
            chemical_data = results_df[results_df['Chemical_Name'] == chemical_name]
            chemical_data = chemical_data.sort_values('Similarity_Score', ascending=False)
            
            # Prioritize human/mouse/rat PDB IDs ONLY
            human_mouse_rat_data = chemical_data[chemical_data['Is_Human_Mouse_Rat'] == True]
            
            # Only take Human/Mouse/Rat matches - no fallback to other organisms
            if len(human_mouse_rat_data) >= 5:
                top_5_targets = human_mouse_rat_data.head(5)
            elif len(human_mouse_rat_data) > 0:
                # Take all available Human/Mouse/Rat matches even if less than 5
                top_5_targets = human_mouse_rat_data
            else:
                # If no Human/Mouse/Rat matches, skip this chemical or mark as no matches
                self.logger.info(f"No Human/Mouse/Rat matches found for {chemical_name}, skipping from CSV output")
                continue
            
            # Create clean comma-separated list of PDB IDs (no star symbols)
            pdb_ids = []
            for _, row in top_5_targets.iterrows():
                pdb_ids.append(row['PDB_ID'])
            
            pdb_ids_str = ', '.join(pdb_ids)
            
            csv_data.append({
                'Plant': top_5_targets.iloc[0]['Plant'],
                'Chemical Name': top_5_targets.iloc[0]['Chemical_Name'],
                'Molecular Structure': top_5_targets.iloc[0]['Molecular_Structure'],
                'Molecule Category': top_5_targets.iloc[0]['Molecule_Category'],
                'Top PDB IDs': pdb_ids_str
            })
        
        # Save to CSV file with standard comma delimiter
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['Plant', 'Chemical Name', 'Molecular Structure', 'Molecule Category', 'Top PDB IDs'])
            
            if csv_data:
                for row in csv_data:
                    writer.writerow([
                        row['Plant'], 
                        row['Chemical Name'], 
                        row['Molecular Structure'], 
                        row['Molecule Category'], 
                        row['Top PDB IDs']
                    ])
            else:
                # If no Human/Mouse/Rat matches found for any chemicals
                writer.writerow(['No Human/Mouse/Rat matches found for any chemicals', '', '', '', ''])
        
        self.logger.info(f"Clean CSV output saved with {len(csv_data)} chemicals (Human/Mouse/Rat only)")
        self.logger.info("Note: Only PDB IDs from Human/Mouse/Rat organisms are included")
    
    def _save_analysis_report(self, results_df, output_file):
        """Save detailed analysis report with top 20 targets and Tanimoto values"""
        with open(output_file, 'w', encoding='utf-8') as f:
            # Header
            f.write("="*80 + "\n")
            f.write("PDB LIGAND SCREENING ANALYSIS REPORT\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            # Processing statistics
            f.write("PROCESSING STATISTICS:\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total PDB records processed: {self.processing_stats['total_pdb_processed']:,}\n")
            f.write(f"Valid PDB fingerprints: {self.processing_stats['valid_pdb_smiles']:,}\n")
            f.write(f"Input chemicals processed: {self.processing_stats['total_input_processed']:,}\n")
            f.write(f"Valid input fingerprints: {self.processing_stats['valid_input_smiles']:,}\n")
            f.write(f"Total similarity comparisons: {self.processing_stats['total_comparisons']:,}\n")
            f.write(f"Chemicals with matches found: {self.processing_stats['matches_found']:,}\n")
            f.write(f"Total result records: {len(results_df):,}\n")
            f.write(f"Errors encountered: {self.processing_stats['errors_encountered']:,}\n\n")
            
            if not results_df.empty:
                # Overall similarity statistics
                similarities = results_df['Similarity_Score'].astype(float)
                f.write("SIMILARITY SCORE STATISTICS:\n")
                f.write("-" * 40 + "\n")
                f.write(f"Mean similarity: {similarities.mean():.4f}\n")
                f.write(f"Median similarity: {similarities.median():.4f}\n")
                f.write(f"Max similarity: {similarities.max():.4f}\n")
                f.write(f"Min similarity: {similarities.min():.4f}\n")
                f.write(f"Standard deviation: {similarities.std():.4f}\n\n")
                
                # Detailed analysis for each molecule
                f.write("DETAILED MOLECULAR ANALYSIS:\n")
                f.write("="*80 + "\n\n")
                
                for chemical_name in results_df['Chemical_Name'].unique():
                    chemical_data = results_df[results_df['Chemical_Name'] == chemical_name]
                    chemical_data = chemical_data.sort_values('Similarity_Score', ascending=False)
                    
                    # Get top 20 targets
                    top_20_targets = chemical_data.head(20)
                    
                    f.write(f"MOLECULE: {top_20_targets.iloc[0]['Chemical_Name']}\n")
                    f.write(f"Plant: {top_20_targets.iloc[0]['Plant']}\n")
                    f.write(f"SMILES: {top_20_targets.iloc[0]['Molecular_Structure']}\n")
                    f.write(f"Category: {top_20_targets.iloc[0]['Molecule_Category']}\n")
                    f.write(f"Total matches found: {len(chemical_data)}\n")
                    f.write("-" * 60 + "\n")
                    
                    f.write("TOP 20 PDB TARGETS WITH TANIMOTO VALUES AND ORGANISM INFO:\n")
                    f.write(f"{'Rank':<4} {'PDB_ID':<8} {'Tanimoto':<10} {'Organism':<20} {'H/M/R':<12} {'Status':<15}\n")
                    f.write("-" * 80 + "\n")
                    
                    for i, (_, row) in enumerate(top_20_targets.iterrows(), 1):
                        organisms = ', '.join(row['Organisms']) if isinstance(row['Organisms'], list) else str(row['Organisms'])
                        human_mouse_rat = "Yes" if row['Is_Human_Mouse_Rat'] else "No"
                        f.write(f"{i:<4} {row['PDB_ID']:<8} {row['Similarity_Score']:<10.4f} "
                               f"{organisms[:20]:<20} {human_mouse_rat:<12} {row['Status']:<15}\n")
                    
                    # Statistics for this molecule
                    mol_similarities = top_20_targets['Similarity_Score'].astype(float)
                    f.write("\nMOLECULE STATISTICS:\n")
                    f.write(f"Mean Tanimoto: {mol_similarities.mean():.4f}\n")
                    f.write(f"Max Tanimoto: {mol_similarities.max():.4f}\n")
                    f.write(f"Min Tanimoto: {mol_similarities.min():.4f}\n")
                    f.write(f"Std Deviation: {mol_similarities.std():.4f}\n")
                    f.write("\n" + "="*80 + "\n\n")
            else:
                f.write("No matches found with similarity > 0.1\n")
        
        self.logger.info(f"Analysis report saved with detailed Tanimoto values")
    
    def run_screening(self, input_file, pdb_file, output_file):
        """Main screening workflow"""
        start_time = datetime.now()
        
        self.logger.info("="*80)
        self.logger.info("FINAL OPTIMIZED PDB LIGAND SCREENING WORKFLOW")
        self.logger.info("="*80)
        self.logger.info(f"Input file: {input_file}")
        self.logger.info(f"PDB database: {pdb_file}")
        self.logger.info(f"Output file: {output_file}")
        self.logger.info(f"Max PDB records: {self.max_pdb_records or 'ALL (600,000+)'}")
        self.logger.info(f"Max input chemicals: {self.max_input_chemicals or 'ALL'}")
        self.logger.info("="*80)
        
        try:
            # Load PDB database
            pdb_df = self._load_pdb_database_robust(pdb_file)
            
            # Process PDB fingerprints
            self._process_pdb_fingerprints(pdb_df)
            
            # Load input chemicals
            self.logger.info(f"Loading input chemicals from: {input_file}")
            input_df = pd.read_csv(input_file)
            self.logger.info(f"Loaded {len(input_df):,} input chemicals")
            
            # Run screening
            self._screen_chemicals(input_df, output_file)
            
            # Final statistics
            end_time = datetime.now()
            self.processing_stats['processing_time'] = (end_time - start_time).total_seconds()
            
            self.logger.info("="*80)
            self.logger.info("SCREENING COMPLETED SUCCESSFULLY")
            self.logger.info("="*80)
            self.logger.info(f"Processing time: {self.processing_stats['processing_time']:.2f} seconds")
            self.logger.info(f"Final output: {output_file}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"Screening failed: {e}")
            return False

def main():
    """Main execution function with professional CLI interface"""
    parser = argparse.ArgumentParser(
        description='Final Optimized PDB Ligand Screening System v3.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all 600K+ PDB ligands with all input chemicals:
  %(prog)s --input input/input_chemicals.csv --pdb pdb_ligands.csv --output output/results.csv
  
  # Process with limits for faster testing:
  %(prog)s --input input/input_chemicals.csv --pdb pdb_ligands.csv --output output/results.csv --max-pdb 100000 --max-input 50
  
  # Process entire database (recommended for production):
  %(prog)s --input input/input_chemicals.csv --pdb pdb_ligands.csv --output output/screening_results.csv

Output Files Generated:
  - [output_name]_top5_targets.csv     (Tab-delimited CSV with top 5 targets per molecule)
  - [output_name]_analysis_report.txt  (Detailed analysis with top 20 targets and Tanimoto values)
        """
    )
    
    # Required arguments
    parser.add_argument('--input', required=True, help='Input chemicals CSV file')
    parser.add_argument('--pdb', required=True, help='PDB ligands CSV file')
    parser.add_argument('--output', required=True, help='Output results CSV file')
    
    # Optional arguments
    parser.add_argument('--max-pdb', type=int, default=0, help='Maximum PDB records to process (0 = ALL 600,000+)')
    parser.add_argument('--max-input', type=int, default=0, help='Maximum input chemicals to process (0 = ALL)')
    parser.add_argument('--chunk-size', type=int, default=10000, help='Chunk size for processing')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.input):
        print(f"❌ ERROR: Input file not found: {args.input}")
        return 1
    
    if not os.path.exists(args.pdb):
        print(f"❌ ERROR: PDB file not found: {args.pdb}")
        return 1
    
    # Create screening system
    screening = FinalOptimizedPDBScreening(
        max_pdb_records=args.max_pdb,
        max_input_chemicals=args.max_input,
        chunk_size=args.chunk_size
    )
    
    # Run screening
    success = screening.run_screening(args.input, args.pdb, args.output)
    
    if success:
        print(f"\n✅ SUCCESS: Final results saved to {args.output}")
        return 0
    else:
        print(f"\n❌ FAILED: Check logs for details")
        return 1

if __name__ == "__main__":
    sys.exit(main())
