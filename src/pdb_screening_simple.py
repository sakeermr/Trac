#!/usr/bin/env python3
"""
PDB Ligand Screening Algorithm for Molecular Similarity Analysis

This module provides a comprehensive molecular similarity screening pipeline
that compares input chemicals against a large database of PDB ligands using
Morgan fingerprints and Tanimoto similarity calculations.

Key Features:
    - Processes 700K+ PDB ligands from the Protein Data Bank
    - Uses RDKit for molecular fingerprint generation
    - Implements chunked processing for large datasets
    - Optimized for GitHub Actions CI/CD environments
    - Robust error handling and logging
    - Memory-efficient processing

Usage:
    python src/pdb_screening_simple.py --input input.csv --pdb pdb_ligands.csv --output results.csv

Requirements:
    - RDKit (conda install -c conda-forge rdkit)
    - pandas >= 2.0.0
    - numpy >= 1.20.0

Author: Standard Seed Corporation
License: See LICENSE file
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

# Suppress RDKit warnings to keep logs clean
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', module='rdkit')

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, DataStructs
    from rdkit import RDLogger
    RDKIT_AVAILABLE = True
    
    # Set up RDKit to be less verbose
    RDLogger.DisableLog('rdApp.*')
    
    # Check if modern fingerprint API is available
    try:
        # Try to create a fingerprint to check if warnings can be suppressed
        from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
        MORGAN_GENERATOR = None
        USE_NEW_API = False
        
    except ImportError:
        # Fall back to old API if new one is not available
        MORGAN_GENERATOR = None
        USE_NEW_API = False
        
except ImportError:
    print("Warning: RDKit not available. Using fallback similarity method.")
    RDKIT_AVAILABLE = False
    MORGAN_GENERATOR = None
    USE_NEW_API = False

# Configure logging to reduce verbose output
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress RDKit logging
logging.getLogger('rdkit').setLevel(logging.ERROR)
logging.getLogger('rdkit.Chem').setLevel(logging.ERROR)
logging.getLogger('rdkit.Chem.rdMolDescriptors').setLevel(logging.ERROR)

def detect_file_encoding(file_path):
    """Detect file encoding by trying different encodings"""
    encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1', 'utf-16', 'utf-32']
    
    for encoding in encodings:
        try:
            with open(file_path, 'r', encoding=encoding) as f:
                # Try to read first few lines
                for i, line in enumerate(f):
                    if i > 100:  # Read first 100 lines
                        break
                return encoding
        except UnicodeDecodeError:
            continue
        except Exception:
            continue
    
    # If all fail, return latin-1 as fallback (it can decode any byte sequence)
    return 'latin-1'

def read_file_in_segments(file_path, encoding='utf-8', segment_size=10000, skip_segments=0):
    """Reads a file in segments, skipping segments with encoding errors"""
    line_buffer = []
    segment_count = 0
    total_lines = 0
    
    try:
        with open(file_path, 'r', encoding=encoding, errors='replace') as f:
            # Skip segments if requested
            if skip_segments > 0:
                logger.info(f"Skipping {skip_segments} segments...")
                for _ in range(skip_segments * segment_size):
                    next(f, None)
                    total_lines += 1
            
            # Read the file segment by segment
            for line in f:
                line_buffer.append(line)
                total_lines += 1
                
                if len(line_buffer) >= segment_size:
                    segment_count += 1
                    yield segment_count, ''.join(line_buffer)
                    line_buffer = []
                    
            # Return any remaining lines
            if line_buffer:
                segment_count += 1
                yield segment_count, ''.join(line_buffer)
                
    except UnicodeDecodeError as e:
        logger.warning(f"Unicode error at line {total_lines}: {e}")
        if line_buffer:
            segment_count += 1
            # Return what we have so far, with errors replaced
            yield segment_count, ''.join(line_buffer)
    
    except Exception as e:
        logger.error(f"Error reading file: {e}")

class SimplePDBScreening:
    def __init__(self, max_pdb_records=1000, max_input_chemicals=100):
        self.pdb_fingerprints = {}
        self.max_pdb_records = max_pdb_records
        self.max_input_chemicals = max_input_chemicals
        
    def validate_smiles(self, smiles):
        """Validate SMILES string with comprehensive error handling"""
        if not smiles or pd.isna(smiles):
            return False
        try:
            smiles_str = str(smiles).strip()
            
            # Basic checks
            if len(smiles_str) < 3:
                return False
            if smiles_str.isspace():
                return False
            
            # Check for common invalid patterns
            invalid_patterns = ['NO_HETEROATOMS', 'WATER', 'HOH', 'CA', 'NA', 'CL']
            if any(pattern in smiles_str.upper() for pattern in invalid_patterns):
                return False
            
            if RDKIT_AVAILABLE:
                # Enhanced RDKit validation with error suppression
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    # Also suppress RDKit internal warnings/errors
                    from rdkit import RDLogger
                    lg = RDLogger.logger()
                    lg.setLevel(RDLogger.CRITICAL)
                    
                    try:
                        mol = Chem.MolFromSmiles(smiles_str)
                        if mol is None:
                            return False
                        
                        # Additional validation: check for reasonable molecular structure
                        if mol.GetNumAtoms() == 0:
                            return False
                        if mol.GetNumAtoms() > 200:  # Skip very large molecules
                            return False
                            
                        # Try to sanitize the molecule to catch valence errors
                        try:
                            Chem.SanitizeMol(mol)
                            return True
                        except:
                            # If sanitization fails, the molecule has structural issues
                            return False
                            
                    except Exception:
                        # Any RDKit parsing error means invalid SMILES
                        return False
                    finally:
                        # Reset logger level
                        lg.setLevel(RDLogger.WARNING)
            else:
                # Simple validation - check if it contains typical SMILES characters
                valid_chars = set('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789()[]{}=#+-@/\\.')
                return len(smiles_str) > 3 and all(c in valid_chars for c in smiles_str)
        except Exception:
            return False
            
    def smiles_to_fingerprint(self, smiles):
        """Convert SMILES to fingerprint using RDKit with comprehensive error handling"""
        try:
            if RDKIT_AVAILABLE:
                # Use RDKit with comprehensive warning suppression
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    # Also suppress RDKit internal warnings
                    from rdkit import RDLogger
                    lg = RDLogger.logger()
                    lg.setLevel(RDLogger.CRITICAL)
                    
                    try:
                        mol = Chem.MolFromSmiles(str(smiles).strip())
                        if mol is None:
                            return None
                        
                        # Try to sanitize the molecule first
                        try:
                            Chem.SanitizeMol(mol)
                        except Exception:
                            # If sanitization fails, try to clean up the molecule
                            try:
                                mol = Chem.MolFromSmiles(str(smiles).strip(), sanitize=False)
                                if mol is None:
                                    return None
                                # Try basic cleanup
                                Chem.FastFindRings(mol)
                            except Exception:
                                return None
                        
                        # Generate fingerprint
                        fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                        
                        # Reset logger level
                        lg.setLevel(RDLogger.WARNING)
                        
                        return fingerprint
                        
                    except Exception:
                        # Any error in fingerprint generation
                        return None
                    finally:
                        # Always reset logger level
                        lg.setLevel(RDLogger.WARNING)
            else:
                # Fallback: simple character-based hash
                return hash(str(smiles).strip()) % (2**20)
        except Exception:
            return None
            
    def calculate_similarity(self, fp1, fp2):
        """Calculate similarity between fingerprints"""
        if fp1 is None or fp2 is None:
            return 0.0
        try:
            if RDKIT_AVAILABLE:
                return DataStructs.TanimotoSimilarity(fp1, fp2)
            else:
                # Simple XOR-based similarity for hash values
                return 1.0 - (abs(fp1 - fp2) / max(fp1, fp2, 1))
        except:
            return 0.0
            
    def load_pdb_ligands(self, pdb_file):
        """Load and process PDB ligands file with robust CSV parsing and error handling"""
        logger.info(f"Loading PDB ligands from {pdb_file}...")
        try:
            # Detect file encoding first
            detected_encoding = detect_file_encoding(pdb_file)
            logger.info(f"Detected file encoding: {detected_encoding}")
            
            # Read file in chunks to handle large files with encoding fallback
            chunk_size = 5000
            processed_count = 0
            
            # Try detected encoding first, then fallback to others
            encodings = [detected_encoding, 'latin-1', 'cp1252', 'iso-8859-1', 'utf-8']
            # Remove duplicates while preserving order
            encodings = list(dict.fromkeys(encodings))
            
            csv_reader = None
            successful_encoding = None
            
            for encoding in encodings:
                try:
                    # Test if we can read the file with this encoding by reading a small sample
                    test_df = pd.read_csv(pdb_file, nrows=100, encoding=encoding, on_bad_lines='skip')
                    logger.info(f"Successfully tested file with {encoding} encoding")
                    
                    # If test successful, create the chunked reader
                    csv_reader = pd.read_csv(
                        pdb_file, 
                        chunksize=chunk_size, 
                        encoding=encoding,
                        quoting=csv.QUOTE_MINIMAL,
                        on_bad_lines='skip',  # Skip malformed lines
                        low_memory=False
                    )
                    successful_encoding = encoding
                    logger.info(f"Successfully opened file with {encoding} encoding")
                    break
                except UnicodeDecodeError as e:
                    logger.warning(f"Unicode error with {encoding} encoding: {e}")
                    continue
                except Exception as e:
                    logger.warning(f"Error with {encoding} encoding: {e}")
                    continue
            
            if csv_reader is None:
                # Final fallback: try reading with the most permissive encoding
                logger.warning("Standard CSV parsing failed, trying final fallback...")
                try:
                    csv_reader = pd.read_csv(
                        pdb_file, 
                        chunksize=chunk_size, 
                        encoding='latin-1',  # latin-1 can decode any byte sequence
                        on_bad_lines='skip',  # Skip bad lines
                        low_memory=False,
                        dtype=str  # Read everything as strings to avoid type issues
                    )
                    successful_encoding = 'latin-1'
                    logger.info("Successfully opened file with final fallback (latin-1)")
                except Exception as e:
                    logger.error(f"All CSV parsing methods failed: {e}")
                    return False
            
            # Process chunks with comprehensive error handling
            chunk_count = 0
            successful_chunks = 0
            
            # Wrap the chunk iteration in a try-catch to handle encoding errors during iteration
            try:
                for chunk_num, chunk in enumerate(csv_reader):
                    chunk_count += 1
                    try:
                        # Additional safety check: convert problematic columns to strings
                        for col in chunk.columns:
                            if chunk[col].dtype == object:
                                chunk[col] = chunk[col].astype(str, errors='ignore')
                        
                        # Handle different column name variations
                        if 'SMILES' not in chunk.columns:
                            smiles_col = None
                            for col in ['smiles', 'Smiles', 'SMILES', 'canonical_smiles']:
                                if col in chunk.columns:
                                    smiles_col = col
                                    break
                            if smiles_col:
                                chunk = chunk.rename(columns={smiles_col: 'SMILES'})
                        
                        if 'PDB_ID' not in chunk.columns:
                            pdb_col = None
                            for col in ['pdb_id', 'PDB_ID', 'pdbid', 'PDB', 'pdb']:
                                if col in chunk.columns:
                                    pdb_col = col
                                    break
                            if pdb_col:
                                chunk = chunk.rename(columns={pdb_col: 'PDB_ID'})
                        
                        # Process each row with error handling
                        for idx, row in chunk.iterrows():
                            try:
                                # Check if we've reached the maximum (0 means no limit)
                                if self.max_pdb_records > 0 and processed_count >= self.max_pdb_records:
                                    break
                                    
                                smiles = row.get('SMILES', '')
                                pdb_id = str(row.get('PDB_ID', '')).strip().upper()
                                heteroatom_code = str(row.get('Heteroatom_Code', '')).strip()
                                
                            except Exception as e:
                                logger.debug(f"Error processing row {idx} in chunk {chunk_num}: {e}")
                                continue
                            
                            # Skip water molecules, ions, and entries with no heteroatoms
                            if heteroatom_code in ['HOH', 'NO_HETEROATOMS', 'CA', 'NA', 'CL', 'MG', 'ZN', 'FE']:
                                continue
                                
                            # Only process entries with valid SMILES and PDB IDs
                            if pdb_id and smiles and pd.notna(smiles) and len(str(smiles).strip()) > 3:
                                if self.validate_smiles(smiles):
                                    fp = self.smiles_to_fingerprint(smiles)
                                    if fp is not None:
                                        # Use PDB_ID + Heteroatom_Code as unique identifier
                                        unique_id = f"{pdb_id}_{heteroatom_code}"
                                        self.pdb_fingerprints[unique_id] = fp
                                        processed_count += 1
                                        
                                        if processed_count % 100 == 0:
                                            logger.info(f"Loaded {processed_count} PDB ligands...")
                                else:
                                    logger.debug(f"Invalid SMILES for {pdb_id}_{heteroatom_code}: {smiles[:50]}...")
                        
                        successful_chunks += 1
                        
                        # Check if we've reached the maximum (0 means no limit)
                        if self.max_pdb_records > 0 and processed_count >= self.max_pdb_records:
                            break
                            
                    except UnicodeDecodeError as e:
                        logger.warning(f"Unicode error in chunk {chunk_num}: {e} - skipping chunk")
                        continue
                    except Exception as e:
                        logger.warning(f"Error processing chunk {chunk_num}: {e} - skipping chunk")
                        continue
                        
            except UnicodeDecodeError as e:
                logger.warning(f"Unicode error during chunk iteration at chunk {chunk_count}: {e}")
                logger.warning(f"Successfully processed {successful_chunks} chunks before error")
                # Don't fail completely - we may have processed some data successfully
            except Exception as e:
                logger.warning(f"Error during chunk iteration at chunk {chunk_count}: {e}")
                logger.warning(f"Successfully processed {successful_chunks} chunks before error")
            
            logger.info(f"Processed {chunk_count} chunks, loaded {processed_count} valid PDB ligands")
            
            # If we got very few results, try alternative processing with smaller chunks
            if processed_count < 100 and chunk_count > 0:
                logger.warning(f"Low yield ({processed_count} ligands from {chunk_count} chunks), trying alternative processing...")
                
                try:
                    # Try with much smaller chunk size to avoid problematic data sections
                    alternative_reader = pd.read_csv(
                        pdb_file, 
                        chunksize=500,  # Much smaller chunks
                        encoding='latin-1',  # Use the most permissive encoding
                        on_bad_lines='skip',
                        low_memory=False,
                        dtype=str,
                        error_bad_lines=False  # Old parameter for compatibility
                    )
                    
                    logger.info("Attempting alternative processing with smaller chunks...")
                    alt_processed = 0
                    
                    for alt_chunk_num, alt_chunk in enumerate(alternative_reader):
                        if alt_chunk_num > 1000:  # Limit to avoid infinite processing
                            break
                            
                        try:
                            # Same processing logic but with more defensive programming
                            for idx, row in alt_chunk.iterrows():
                                try:
                                    smiles = str(row.get('SMILES', '')).strip()
                                    pdb_id = str(row.get('PDB_ID', '')).strip().upper()
                                    heteroatom_code = str(row.get('Heteroatom_Code', '')).strip()
                                    
                                    if (pdb_id and smiles and len(smiles) > 3 and 
                                        heteroatom_code not in ['HOH', 'NO_HETEROATOMS', 'CA', 'NA', 'CL', 'MG', 'ZN', 'FE']):
                                        
                                        if self.validate_smiles(smiles):
                                            fp = self.smiles_to_fingerprint(smiles)
                                            if fp is not None:
                                                unique_id = f"{pdb_id}_{heteroatom_code}"
                                                if unique_id not in self.pdb_fingerprints:
                                                    self.pdb_fingerprints[unique_id] = fp
                                                    alt_processed += 1
                                                    
                                                    if alt_processed % 100 == 0:
                                                        logger.info(f"Alternative processing: {alt_processed} additional ligands...")
                                                        
                                except Exception:
                                    continue
                                    
                        except Exception as e:
                            logger.debug(f"Error in alternative chunk {alt_chunk_num}: {e}")
                            continue
                    
                    if alt_processed > 0:
                        logger.info(f"Alternative processing added {alt_processed} additional ligands")
                        processed_count += alt_processed
                        
                except Exception as e:
                    logger.warning(f"Alternative processing failed: {e}")
            
            return processed_count > 0
            
        except UnicodeDecodeError as e:
            logger.warning(f"Final encoding error loading PDB file: {e}")
            logger.warning("The file contains corrupted or mixed encoding data")
            
            # If we have some data already loaded, return partial success
            if processed_count > 0:
                logger.info(f"Partial success: loaded {processed_count} valid PDB ligands before encoding error")
                return True
            else:
                logger.error("No ligands could be loaded due to encoding issues")
                logger.error("Try converting the file to UTF-8 encoding or use a different file format")
                return False
            return False
        except Exception as e:
            logger.error(f"General error loading PDB file: {e}")
            
            # If we have some data already loaded, return partial success
            if processed_count > 0:
                logger.info(f"Partial success: loaded {processed_count} valid PDB ligands before error")
                return True
            else:
                logger.error("Failed to load PDB ligands")
                return False
            return False
            
    def _process_file_segments(self, file_path, start_segment=0, max_segments=None):
        """Process file in segments to bypass encoding errors"""
        valid_entries = 0
        total_segments = 0
        segment_size = 10000
        processed_count = 0
        
        logger.info(f"Starting segment-based processing from segment {start_segment}...")
        
        # Process each segment
        for segment_num, segment_data in read_file_in_segments(file_path, encoding='latin-1', 
                                                            segment_size=segment_size, skip_segments=start_segment):
            if max_segments and segment_num > start_segment + max_segments:
                break
                
            total_segments += 1
            
            try:
                # Parse segment as CSV
                import io
                segment_df = pd.read_csv(io.StringIO(segment_data), on_bad_lines='skip')
                
                # Handle different column name variations
                if 'SMILES' not in segment_df.columns:
                    smiles_col = None
                    for col in ['smiles', 'Smiles', 'SMILES', 'canonical_smiles']:
                        if col in segment_df.columns:
                            smiles_col = col
                            break
                    if smiles_col:
                        segment_df = segment_df.rename(columns={smiles_col: 'SMILES'})
                
                if 'PDB_ID' not in segment_df.columns:
                    pdb_col = None
                    for col in ['pdb_id', 'PDB_ID', 'pdbid', 'PDB', 'pdb']:
                        if col in segment_df.columns:
                            pdb_col = col
                            break
                    if pdb_col:
                        segment_df = segment_df.rename(columns={pdb_col: 'PDB_ID'})
                
                # Process each row with error handling
                valid_in_segment = 0
                for idx, row in segment_df.iterrows():
                    try:
                        # Check if we've reached the maximum (0 means no limit)
                        if self.max_pdb_records > 0 and processed_count >= self.max_pdb_records:
                            break
                            
                        smiles = row.get('SMILES', '')
                        pdb_id = str(row.get('PDB_ID', '')).strip().upper()
                        heteroatom_code = str(row.get('Heteroatom_Code', '')).strip()
                        
                    except Exception as e:
                        logger.debug(f"Error processing row {idx} in segment {segment_num}: {e}")
                        continue
                    
                    # Skip water molecules, ions, and entries with no heteroatoms
                    if heteroatom_code in ['HOH', 'NO_HETEROATOMS', 'CA', 'NA', 'CL', 'MG', 'ZN', 'FE']:
                        continue
                        
                    # Only process entries with valid SMILES and PDB IDs
                    if pdb_id and smiles and pd.notna(smiles) and len(str(smiles).strip()) > 3:
                        if self.validate_smiles(smiles):
                            fp = self.smiles_to_fingerprint(smiles)
                            if fp is not None:
                                # Use PDB_ID + Heteroatom_Code as unique identifier
                                unique_id = f"{pdb_id}_{heteroatom_code}"
                                self.pdb_fingerprints[unique_id] = fp
                                processed_count += 1
                                valid_in_segment += 1
                                
                                if processed_count % 1000 == 0:
                                    logger.info(f"Loaded {processed_count} PDB ligands (segment {segment_num})...")
                        else:
                            logger.debug(f"Invalid SMILES for {pdb_id}_{heteroatom_code}: {smiles[:50]}...")
                
                valid_entries += valid_in_segment
                
                # Log progress periodically
                if segment_num % 5 == 0 or valid_in_segment > 0:
                    logger.info(f"Processed segment {segment_num}: {valid_in_segment} valid entries ({valid_entries} total)")
                    
            except Exception as e:
                logger.warning(f"Error processing segment {segment_num}: {e}")
                continue
        
        logger.info(f"Segment processing complete: {total_segments} segments, {valid_entries} valid entries")
        return valid_entries
        
    def load_pdb_ligands_with_segments(self, pdb_file, start_segment=0, max_segments=None):
        """Load PDB ligands using segment-based approach to bypass encoding errors"""
        logger.info(f"Loading PDB ligands with segment-based processing from {pdb_file}...")
        
        # Start with empty fingerprints if loading from scratch
        if start_segment == 0:
            self.pdb_fingerprints = {}
        
        try:
            # Process file in segments
            valid_entries = self._process_file_segments(pdb_file, start_segment, max_segments)
            
            if valid_entries > 0:
                logger.info(f"Successfully loaded {valid_entries} PDB ligands using segment-based processing")
                logger.info(f"Total fingerprints in memory: {len(self.pdb_fingerprints)}")
                return True
            else:
                logger.error("Failed to load any valid PDB ligands with segment-based processing")
                return False
        
        except Exception as e:
            logger.error(f"Error in segment-based PDB loading: {e}")
            # If we have some data already loaded, return partial success
            if len(self.pdb_fingerprints) > 0:
                logger.info(f"Partial success: {len(self.pdb_fingerprints)} PDB ligands in memory")
                return True
            else:
                return False

    def find_similar_ligands(self, target_smiles, threshold=0.5, top_n=5):
        """Find similar ligands in PDB database"""
        target_fp = self.smiles_to_fingerprint(target_smiles)
        if target_fp is None:
            return []
            
        similarities = []
        for pdb_id, pdb_fp in self.pdb_fingerprints.items():
            sim = self.calculate_similarity(target_fp, pdb_fp)
            if sim >= threshold:
                similarities.append((pdb_id, sim))
                
        similarities.sort(key=lambda x: x[1], reverse=True)
        return similarities[:top_n]
        
    def process_chemicals(self, input_file, output_file):
        """Process input chemicals and find PDB matches with robust encoding handling"""
        logger.info(f"Processing chemicals from {input_file}")
        
        try:
            # Read input file with encoding fallback
            df = None
            encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(input_file, encoding=encoding)
                    logger.info(f"Successfully opened input file with {encoding} encoding")
                    break
                except UnicodeDecodeError:
                    logger.warning(f"Failed to read input file with {encoding} encoding, trying next...")
                    continue
                except Exception as e:
                    logger.warning(f"Error with {encoding} encoding: {e}")
                    continue
            
            if df is None:
                logger.error("Failed to read input file with any encoding")
                return
            
            # Standardize column names
            df.columns = df.columns.str.lower().str.replace(' ', '.')
            
            # Limit input size if specified (0 means no limit)
            if self.max_input_chemicals > 0 and len(df) > self.max_input_chemicals:
                logger.warning(f"Limiting input to {self.max_input_chemicals} chemicals")
                df = df.head(self.max_input_chemicals)
            else:
                logger.info(f"Processing all {len(df)} input chemicals (no limit set)")
            
            # Initialize result columns (only keeping pdb_targets)
            df['pdb_targets'] = ''
            # Temporary column for processing status (will be removed before saving)
            df['_temp_processing_status'] = 'pending'
            
            processed = 0
            successful = 0
            
            for idx, row in df.iterrows():
                try:
                    # Get SMILES from various possible column names
                    smiles = None
                    for col in ['molecular.structure', 'smiles', 'structure', 'mol_structure']:
                        if col in df.columns and pd.notna(row.get(col)):
                            smiles = str(row[col]).strip()
                            break
                    
                    if not smiles or not self.validate_smiles(smiles):
                        df.at[idx, '_temp_processing_status'] = 'invalid_smiles'
                        processed += 1
                        continue
                    
                    # Find similar ligands with lower threshold for better results
                    similar_ligands = self.find_similar_ligands(smiles, threshold=0.2, top_n=5)
                    
                    if similar_ligands:
                        # Extract only PDB IDs (before underscore) and remove duplicates
                        pdb_ids = []
                        for lig_id, score in similar_ligands:
                            # Extract PDB ID (before underscore)
                            pdb_id = lig_id.split('_')[0]
                            if pdb_id not in pdb_ids:  # Avoid duplicates
                                pdb_ids.append(pdb_id)
                        
                        # Join with comma and space (cleaner format)
                        df.at[idx, 'pdb_targets'] = ', '.join(pdb_ids)
                        df.at[idx, '_temp_processing_status'] = 'success'
                        successful += 1
                    else:
                        df.at[idx, '_temp_processing_status'] = 'no_matches'
                    
                    processed += 1
                    
                    if processed % 10 == 0:
                        logger.info(f"Processed {processed}/{len(df)} chemicals")
                        
                except Exception as e:
                    logger.error(f"Error processing row {idx}: {e}")
                    df.at[idx, '_temp_processing_status'] = f'error: {str(e)[:50]}'
                    processed += 1
            
            # Remove temporary processing status column before saving
            if '_temp_processing_status' in df.columns:
                df.drop(columns=['_temp_processing_status'], inplace=True)
            
            # Save results (only original columns + pdb_targets)
            df.to_csv(output_file, index=False, encoding='utf-8')
            logger.info(f"Results saved to {output_file}")
            logger.info(f"Total processed: {processed}, Successful: {successful}")
            
            return processed, successful, output_file
            
        except UnicodeDecodeError as e:
            logger.error(f"Encoding error processing chemicals: {e}")
            logger.error("Try converting the input file to UTF-8 encoding")
            return 0, 0, None
        except Exception as e:
            logger.error(f"Error processing chemicals: {e}")
            return 0, 0, None

def main():
    parser = argparse.ArgumentParser(description='PDB Ligand Screening')
    parser.add_argument('--input', required=True, help='Input CSV file')
    parser.add_argument('--pdb', required=True, help='PDB ligands CSV file')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--max-pdb', type=int, default=1000, help='Max PDB records to process')
    parser.add_argument('--max-input', type=int, default=100, help='Max input chemicals to process')
    
    args = parser.parse_args()
    
    # Initialize screening
    screening = SimplePDBScreening(
        max_pdb_records=args.max_pdb,
        max_input_chemicals=args.max_input
    )
    
    # Load PDB ligands
    if not screening.load_pdb_ligands(args.pdb):
        logger.error("Failed to load PDB ligands")
        sys.exit(1)
    
    # Process chemicals
    processed, successful, output_file = screening.process_chemicals(args.input, args.output)
    
    if output_file:
        print(f"SUCCESS: Processed {processed} chemicals, {successful} successful matches")
        print(f"Output saved to: {output_file}")
        
        # Output for GitHub Actions (handle Windows/Unix differences)
        github_output = os.environ.get('GITHUB_OUTPUT', None)
        if github_output:
            try:
                with open(github_output, 'a') as f:
                    f.write(f"total_processed={processed}\n")
                    f.write(f"successful={successful}\n")
                    f.write(f"output_file={output_file}\n")
            except:
                # If GitHub Actions output file is not accessible, skip silently
                pass
    else:
        logger.error("Processing failed")
        sys.exit(1)

if __name__ == "__main__":
    main()
