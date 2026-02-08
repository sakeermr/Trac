#!/usr/bin/env python3
"""
Input Validation Module for PDB Ligand Screening Pipeline

This module provides comprehensive validation of input CSV files to ensure
data quality before molecular similarity analysis. It validates chemical
names, SMILES structures, and file formats.

Key Features:
    - CSV format validation
    - SMILES string validation
    - Chemical name verification
    - File existence checks
    - Data type validation
    - Comprehensive error reporting

Usage:
    python src/validate_input.py --input input/input_chemicals.csv

Requirements:
    - pandas >= 2.0.0
    - Standard Python libraries

Author: Standard Seed Corporation
License: See LICENSE file
"""

import pandas as pd
import sys
import os
from typing import Dict, List, Tuple

def validate_smiles(smiles: str) -> bool:
    """
    Validate SMILES string using basic checks.
    For full validation, RDKit would be needed.
    """
    if not smiles or pd.isna(smiles):
        return False
    
    smiles = str(smiles).strip()
    if not smiles:
        return False
    
    # Basic SMILES validation (without RDKit)
    # Check for obviously invalid characters
    invalid_chars = ['@', '#', '$', '%', '^', '&', '*', '!', '?', '<', '>', '|', '\\', '/', '`', '~']
    for char in invalid_chars:
        if char in smiles:
            return False
    
    # Must contain at least one carbon or nitrogen
    if not any(char in smiles.upper() for char in ['C', 'N', 'O', 'S', 'P']):
        return False
    
    return True

def validate_csv_format(file_path: str) -> Dict:
    """
    Validate CSV file format and content.
    
    Returns:
        Dict with validation results
    """
    results = {
        'valid': False,
        'file_exists': False,
        'readable': False,
        'has_header': False,
        'has_required_columns': False,
        'total_rows': 0,
        'valid_smiles_count': 0,
        'errors': [],
        'warnings': [],
        'column_mapping': {},
        'sample_data': []
    }
    
    # Check if file exists
    if not os.path.exists(file_path):
        results['errors'].append(f"File not found: {file_path}")
        return results
    
    results['file_exists'] = True
    
    # Try to read the file
    try:
        df = pd.read_csv(file_path)
        results['readable'] = True
        results['total_rows'] = len(df)
    except Exception as e:
        results['errors'].append(f"Error reading CSV file: {e}")
        return results
    
    # Check if file is empty
    if len(df) == 0:
        results['errors'].append("CSV file is empty")
        return results
    
    # Get column names
    columns = [col.strip() for col in df.columns]
    results['has_header'] = len(columns) > 0
    
    # Define required columns and possible variations
    required_columns = {
        'Chemical Name': ['chemical name', 'chemical_name', 'name', 'compound_name', 'compound name'],
        'Molecular Structure': ['molecular structure', 'molecular_structure', 'smiles', 'structure'],
        'Molecule Category': ['molecule category', 'molecule_category', 'category', 'type', 'class']
    }
    
    # Find column mappings
    column_mapping = {}
    for required_col, variations in required_columns.items():
        found = False
        for col in columns:
            if col.lower() in variations:
                column_mapping[required_col] = col
                found = True
                break
        
        if not found:
            results['errors'].append(f"Required column '{required_col}' not found. Looking for one of: {variations}")
    
    results['column_mapping'] = column_mapping
    results['has_required_columns'] = len(column_mapping) == 3
    
    if not results['has_required_columns']:
        results['errors'].append(f"Available columns: {columns}")
        return results
    
    # Rename columns for easier processing
    df_renamed = df.rename(columns={v: k for k, v in column_mapping.items()})
    
    # Validate data content
    valid_smiles_count = 0
    invalid_rows = []
    
    for idx, row in df_renamed.iterrows():
        chemical_name = row.get('Chemical Name', '')
        molecular_structure = row.get('Molecular Structure', '')
        molecule_category = row.get('Molecule Category', '')
        
        # Check for missing data
        if pd.isna(chemical_name) or str(chemical_name).strip() == '':
            invalid_rows.append(f"Row {idx + 2}: Missing chemical name")
        
        if pd.isna(molecular_structure) or str(molecular_structure).strip() == '':
            invalid_rows.append(f"Row {idx + 2}: Missing molecular structure")
        else:
            # Validate SMILES
            if validate_smiles(molecular_structure):
                valid_smiles_count += 1
            else:
                invalid_rows.append(f"Row {idx + 2}: Invalid SMILES '{molecular_structure}'")
        
        if pd.isna(molecule_category) or str(molecule_category).strip() == '':
            invalid_rows.append(f"Row {idx + 2}: Missing molecule category")
    
    results['valid_smiles_count'] = valid_smiles_count
    
    # Add warnings for invalid rows
    if invalid_rows:
        results['warnings'].extend(invalid_rows[:10])  # Show first 10 issues
        if len(invalid_rows) > 10:
            results['warnings'].append(f"... and {len(invalid_rows) - 10} more issues")
    
    # Create sample data preview
    if len(df_renamed) > 0:
        sample_size = min(3, len(df_renamed))
        for idx in range(sample_size):
            row = df_renamed.iloc[idx]
            results['sample_data'].append({
                'Chemical Name': row.get('Chemical Name', ''),
                'Molecular Structure': row.get('Molecular Structure', ''),
                'Molecule Category': row.get('Molecule Category', '')
            })
    
    # Overall validation
    results['valid'] = (
        results['file_exists'] and 
        results['readable'] and 
        results['has_required_columns'] and 
        results['valid_smiles_count'] > 0
    )
    
    return results

def print_validation_report(results: Dict) -> None:
    """Print a formatted validation report."""
    print("🔍 CSV Validation Report")
    print("=" * 50)
    
    # File status
    if results['file_exists']:
        print("✅ File exists")
    else:
        print("❌ File not found")
        return
    
    if results['readable']:
        print("✅ File is readable")
        print(f"   Total rows: {results['total_rows']}")
    else:
        print("❌ File is not readable")
        return
    
    # Column validation
    if results['has_required_columns']:
        print("✅ Required columns found")
        for standard_name, actual_name in results['column_mapping'].items():
            print(f"   {standard_name} → {actual_name}")
    else:
        print("❌ Missing required columns")
    
    # SMILES validation
    if results['valid_smiles_count'] > 0:
        print(f"✅ Valid SMILES found: {results['valid_smiles_count']}/{results['total_rows']}")
        success_rate = (results['valid_smiles_count'] / results['total_rows']) * 100
        print(f"   Success rate: {success_rate:.1f}%")
    else:
        print("❌ No valid SMILES found")
    
    # Show sample data
    if results['sample_data']:
        print("\n📋 Sample Data (first 3 rows):")
        for i, row in enumerate(results['sample_data'], 1):
            print(f"Row {i}:")
            print(f"  Name: {row['Chemical Name']}")
            print(f"  SMILES: {row['Molecular Structure']}")
            print(f"  Category: {row['Molecule Category']}")
            print()
    
    # Show errors
    if results['errors']:
        print("❌ Errors:")
        for error in results['errors']:
            print(f"   {error}")
        print()
    
    # Show warnings
    if results['warnings']:
        print("⚠️ Warnings:")
        for warning in results['warnings']:
            print(f"   {warning}")
        print()
    
    # Overall status
    if results['valid']:
        print("🎉 Overall Status: VALID - Ready for processing!")
    else:
        print("❌ Overall Status: INVALID - Please fix issues before processing")

def main():
    """Main validation function."""
    if len(sys.argv) != 2:
        print("Usage: python validate_input.py <input_csv_file>")
        print("Example: python validate_input.py input/input_chemicals.csv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    print(f"Validating: {input_file}")
    print()
    
    results = validate_csv_format(input_file)
    print_validation_report(results)
    
    # Exit with appropriate code
    if results['valid']:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == "__main__":
    main()
