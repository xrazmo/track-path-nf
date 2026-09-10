import pandas as pd
import numpy as np
import os
import io
import re
import concurrent.futures
from collections import defaultdict
from tqdm import tqdm
import time

def parse_vcf_file(file_path):
    """
    Optimized function to parse a VCF file and convert it to a pandas DataFrame
    
    Args:
        file_path (str): Path to the VCF file
        
    Returns:
        pandas.DataFrame: DataFrame containing variant information
    """
    # Extract sample ID from filename
    sample_id = os.path.basename(file_path).split('.')[0]
    
    # Initialize data containers
    header = None
    data_rows = []
    
    # Process file with buffer reading for memory efficiency
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('##'):
                continue  # Skip metadata lines
            elif line.startswith('#'):
                # Extract column headers
                header = line[1:].strip().split('\t')
            else:
                # Store data lines directly as lists
                data_rows.append(line.strip().split('\t'))
                
    
    # Validate header
    if not header:
        raise ValueError(f"No column headers found in VCF file: {file_path}")
    
    # Create DataFrame all at once (faster than appending)
    df = pd.DataFrame(data_rows, columns=header)
    
    # Process INFO field efficiently using vectorized operations
    df = process_info_field(df)

    # Process ANN field if available
    if 'INFO_ANN' in df.columns:
        process_ann_field(df)
    
    # Process FORMAT fields efficiently for all samples at once
    if len(header) > 8:  # There are sample columns
        process_format_fields(df, header[9:])
    
    # Add sample ID column
    df['sample_id'] = sample_id
    # Rename column with the value of sample_id to 'samples_quality'
    df.rename(columns={sample_id: 'FORMAT_DATA'}, inplace=True)
     # Drop all columns except the specified ones
    keep_columns = [
        'CHROM', 'sample_id', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FORMAT', 'FORMAT_DATA',
        'INFO_AB', 'INFO_AO', 'INFO_DP', 'INFO_QA', 'INFO_QR', 'INFO_RO', 'INFO_TYPE', 'INFO_LOF',
        'ANN_Allele', 'ANN_Annotation', 'ANN_Annotation_Impact', 'ANN_Gene_Name', 'ANN_Gene_ID',
         'ANN_Feature_ID', 'ANN_HGVS.c','ANN_HGVS.p', 'ANN_cDNA_pos', 'ANN_CDS_pos', 'ANN_AA_pos'
    ]
    # Only keep columns that exist in the DataFrame
    cols_to_keep = [col for col in keep_columns if col in df.columns]
    df.drop(columns=[col for col in df.columns if col not in cols_to_keep], inplace=True)
   
    return df

def process_info_field(df):
    """
    Process INFO field efficiently using vectorized operations
    
    Args:
        df (pandas.DataFrame): DataFrame containing VCF data with INFO column
    """
    
    info_dicts = []
    for idx in range(len(df)):  # Use range to ensure all indices are processed
    # Initialize a dictionary for this row with its index
        row_dict = {'idx': idx}
        
        # Get the INFO value safely
        if idx < len(df) and 'INFO' in df.columns:
            info_str = df['INFO'].iloc[idx]
            
            if pd.notna(info_str):  # Skip NaN values, but keep the row
                try:
                    for field in str(info_str).split(';'):
                        if '=' in field:
                            key, value = field.split('=', 1)
                            row_dict[key] = value
                        
                except Exception as e:
                    print(f"Error processing INFO field at index {idx}: {info_str}, Error: {e}")
        
        # Add the dictionary to our list regardless of whether we found info
        info_dicts.append(row_dict)

    # Create a new DataFrame from the list of dictionaries
    info_df = pd.DataFrame(info_dicts)
    # Add INFO_ prefix to all columns except 'idx'
    info_df = info_df.rename(columns={col: f'INFO_{col}' if col != 'idx' else col for col in info_df.columns})
   
    # Merge the new DataFrame with the original one based on the index
    # Using outer join to ensure no rows are lost
    mdf = pd.merge(
        df, 
        info_df,
        left_index=True, 
        right_on='idx', 
        how='outer'  # Using outer join to ensure all rows from both DataFrames are included
    )
    return mdf


def process_ann_field(df):
    """
    Process ANN field efficiently
    
    Args:
        df (pandas.DataFrame): DataFrame containing VCF data with INFO_ANN column
    """
    # Check if INFO_ANN column exists
    if 'INFO_ANN' not in df.columns:
        return
        
    ann_columns = ['Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID', 
                  'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c', 
                  'HGVS.p', 'cDNA_pos', 'CDS_pos', 'AA_pos', 'Distance', 'ERRORS_WARNINGS_INFO']
    
    # Initialize arrays for each annotation column
    ann_data = {col: np.full(len(df), None, dtype=object) for col in ann_columns}
    
    # Process annotations for rows that have them
    mask = df['INFO_ANN'].notna()
    
    if not mask.any():
        return
    # Process only rows with annotations
    for idx in df.index[mask]:
        try:
            ann_value = df.loc[idx, 'INFO_ANN']
            # Take the first annotation if multiple exist
            first_ann = str(ann_value).split(',')[0]
            ann_parts = first_ann.split('|')
            # Add annotation parts to their respective columns
            for i, col in enumerate(ann_columns):
                if i < len(ann_parts):
                    ann_data[col][idx] = ann_parts[i]
        except Exception as e:
            print(f"Error processing ANN field at index {idx}: {ann_value}, Error: {e}")
            continue
    
    # Add all annotation columns at once
    for col, values in ann_data.items():
        df[f'ANN_{col}'] = values

def process_format_fields(df, sample_columns):
    """
    Process FORMAT fields efficiently for all samples
    
    Args:
        df (pandas.DataFrame): DataFrame containing VCF data
        sample_columns (list): List of sample column names
    """
    # Skip if FORMAT column is empty
    if df.empty or not pd.notna(df['FORMAT'].iloc[0]):
        return
   
    # Check if FORMAT is consistent across all rows
    if df['FORMAT'].nunique() == 1:
        # Fast path: all rows have the same FORMAT
        format_keys = str(df['FORMAT'].iloc[0]).split(':')
        
        for sample in sample_columns:
            # Pre-allocate arrays for each format field
            format_data = {key: np.full(len(df), None, dtype=object) for key in format_keys}
            
            # Process all rows for this sample at once
            for idx, format_value in enumerate(df[sample]):
                if pd.notna(format_value):
                    try:
                        # Ensure format_value is a string before splitting
                        format_values = str(format_value).split(':')
                        for i, key in enumerate(format_keys):
                            if i < len(format_values):
                                format_data[key][idx] = format_values[i]
                    except Exception as e:
                        # Skip problematic values but log for debugging
                        print(f"Error processing format value at index {idx}: {format_value}, Error: {e}")
                        continue
            
            # Add all format columns for this sample at once
            for key, values in format_data.items():
                df[f'{sample}_{key}'] = values
    else:
        # Slow path: FORMAT varies by row
        for sample in sample_columns:
            # Get all unique format keys across all rows
            all_format_keys = set()
            for fmt in df['FORMAT'].unique():
                if pd.notna(fmt):  # Check for NaN values
                    all_format_keys.update(str(fmt).split(':'))
            
            # Pre-allocate arrays for each format field
            format_data = {key: np.full(len(df), None, dtype=object) for key in all_format_keys}
            
            # Process each row with its specific format
            for idx, row in df.iterrows():
                if pd.notna(row['FORMAT']) and pd.notna(row[sample]):
                    try:
                        format_keys = str(row['FORMAT']).split(':')
                        format_values = str(row[sample]).split(':')
                        
                        for i, key in enumerate(format_keys):
                            if i < len(format_values):
                                format_data[key][idx] = format_values[i]
                    except Exception as e:
                        print(f"Error processing row {idx}: FORMAT={row['FORMAT']}, {sample}={row[sample]}, Error: {e}")
                        continue
            
            # Add all format columns for this sample at once
            for key, values in format_data.items():
                df[f'{sample}_{key}'] = values

def process_vcf_files_parallel(vcf_files, max_workers=None):
    """
    Process multiple VCF files in parallel and return a list of DataFrames
    
    Args:
        vcf_files (list): List of paths to VCF files
        max_workers (int, optional): Maximum number of worker processes. 
                                    If None, uses CPU count.
    
    Returns:
        list: List of pandas DataFrames, each containing data from one VCF file
    """
    dataframes = []
    
    # Define a wrapper function to use with tqdm
    def process_file(file_path):
        try:
            return parse_vcf_file(file_path)
        except Exception as e:
            print(f"Error processing {file_path}: {str(e)}")
            return None
    
    # Use ThreadPoolExecutor for I/O bound tasks
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks and wrap with tqdm for progress bar
        futures = {executor.submit(process_file, file): file for file in vcf_files}
        
        # Process results as they complete with progress bar
        for future in tqdm(concurrent.futures.as_completed(futures), 
                          total=len(futures), 
                          desc="Processing VCF files"):
            file = futures[future]
            try:
                df = future.result()
                if df is not None:
                    dataframes.append(df)
            except Exception as e:
                print(f"Error processing {file}: {str(e)}")
    
    return dataframes

def save_dataframes(dataframes, output_dir, prefix="vcf_data"):
    """
    Save a list of DataFrames to CSV files
    
    Args:
        dataframes (list): List of pandas DataFrames
        output_dir (str): Directory to save CSV files
        prefix (str): Prefix for CSV filenames
    """
    os.makedirs(output_dir, exist_ok=True)
    
    for i, df in enumerate(dataframes):
        output_file = os.path.join(output_dir, f"{prefix}_{i}.csv")
        df.to_csv(output_file, index=False)
        print(f"Saved DataFrame {i} to {output_file}")

def save_combined_dataframe(dataframes, output_file):
    """
    Combine DataFrames and save to a single CSV file
    
    Args:
        dataframes (list): List of pandas DataFrames
        output_file (str): Path to the output CSV file
    """
    if not dataframes:
        print("No DataFrames to combine")
        return
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Combine DataFrames
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    # Save combined DataFrame
    combined_df.to_csv(output_file, index=False)
    print(f"Saved combined DataFrame to {output_file}")
    print(f"Combined DataFrame shape: {combined_df.shape}")

# For command line usage
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Process VCF files in parallel')
    parser.add_argument('--input_dir', required=True, help='Directory containing VCF files')
    parser.add_argument('--pattern', default='*.vcf', help='File pattern to match VCF files (default: *.vcf)')
    parser.add_argument('--output_dir', default='output', help='Directory to save output files')
    parser.add_argument('--combine', action='store_true', help='Combine all DataFrames into one CSV file')
    parser.add_argument('--workers', type=int, default=None, help='Maximum number of worker processes')
    
    args = parser.parse_args()
    
    # Get list of VCF files
    import glob
    vcf_files = glob.glob(os.path.join(args.input_dir, args.pattern))
    
    if not vcf_files:
        print(f"No VCF files found matching pattern '{args.pattern}' in '{args.input_dir}'")
        exit(1)
    
    print(f"Found {len(vcf_files)} VCF files to process")
    
    # Process VCF files in parallel
    start_time = time.time()
    dataframes = process_vcf_files_parallel(vcf_files, max_workers=args.workers)
    end_time = time.time()
    
    print(f"Processed {len(dataframes)} VCF files in {end_time - start_time:.2f} seconds")
    
    # Save DataFrames
    if args.combine:
        output_file = os.path.join(args.output_dir, "combined_vcf_data.csv")
        save_combined_dataframe(dataframes, output_file)
    else:
        save_dataframes(dataframes, args.output_dir)
    
    print("Done!")