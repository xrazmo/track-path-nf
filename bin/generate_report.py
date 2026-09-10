import argparse
import copy
import csv
import os
from glob import glob
import shutil
import zipfile
import pandas as pd
import re
import json

def extract_mlst_data(input_dir):

    """
    Extract MLST (Multi-Locus Sequence Typing) data from a specified input directory.
    This function processes MLST data files generated using the `tseemann/mlst` pipeline.
    It assumes the input directory contains a subdirectory named `mlst` with files
    in the format `[sample-id].mlst.tsv`. Each TSV file is expected to contain information
    about the MLST scheme, sequence type (ST), and allele data.
    Args:
        input_dir (str): The path to the input directory containing the `mlst` folder.
    Returns:
        dict: A dictionary where each key is a sample ID and the value is another dictionary
            containing the following keys:
            - "scheme" (str): The MLST scheme name.
            - "st" (str): The sequence type (ST).
            - "loci" (list): A list of loci names.
            - "alleles" (list): A list of alleles corresponding to the loci. Alleles are
                integers if they can be converted, otherwise they remain as strings.
    Notes:
        - If the `mlst` directory does not exist or no TSV files are found, the function
        prints an error message and returns an empty list.
        - Allele data is extracted using a regex pattern to match the format `locus(allele)`.
        - Question marks in allele values are removed, and the function attempts to convert
        allele values to integers where possible.
    """
   
    mlst_dir = os.path.join(input_dir, "mlst")
    if not os.path.exists(mlst_dir):
        print(f"Error: MLST directory '{mlst_dir}' does not exist.")
        return {}

    mlst_tsv_files = glob(os.path.join(mlst_dir, "*.tsv"))
    if not mlst_tsv_files:
        print(f"Error: No MLST TSV files found in '{mlst_dir}'.")
        return {}

    mlst_data = {}

    for filename in mlst_tsv_files:
        sample_id = os.path.basename(filename).split('.')[0]
        mlst_data[sample_id] = {
            "scheme": "-",
            "st": "-",
            "loci": [],
            "alleles": []
        }

        with open(filename, 'r') as file:
            for line in file:
                parts = line.strip().split("\t")

                if len(parts) > 1:
                    mlst_data[sample_id]["scheme"] = parts[1]

                if len(parts) > 2:
                    mlst_data[sample_id]["st"] = parts[2]

                if len(parts) > 3:
                    for allele_info in parts[3:]:
                        match = re.match(r'([^(]+)\(([^)]+)\)', allele_info)
                        if match:
                            locus = match.group(1)
                            allele = match.group(2)

                            mlst_data[sample_id]["loci"].append(locus)
                            try:
                                mlst_data[sample_id]["alleles"].append(int(allele))
                            except ValueError:
                                mlst_data[sample_id]["alleles"].append(allele)

    return mlst_data

def extract_bracken_data(input_dir):
    """
    Extracts and processes Bracken data from a specified input directory.
    This function reads `.bracken.tsv` files from the `kraken2_bracken` subdirectory
    within the given `input_dir`. It processes each file to extract the top taxonomy
    entry based on the highest fraction of total reads and returns the data in a
    structured dictionary format.
    Args:
        input_dir (str): Path to the input directory containing the `kraken2_bracken` folder.
    Returns:
        dict: A dictionary where keys are sample IDs (derived from file names) and values
              are dictionaries containing the top taxonomy entry with the following fields:
              - 'name': Taxonomy name.
              - 'taxid': Taxonomy ID.
              - 'lvl': Taxonomy level.
              - 'ar': Added reads.
              - 'nr': New estimated reads.
              - 'kar': Kraken-assigned reads.
              - 'frac': Fraction of total reads (rounded to 2 decimal places).
    Notes:
        - If the `kraken2_bracken` directory does not exist, an error message is printed,
          and an empty dictionary is returned.
        - If a `.bracken.tsv` file is empty, a warning is printed, and the file is skipped.
        - The function assumes that the `.bracken.tsv` files have a specific column structure
          that matches the `column_mapping` dictionary.
    Example:
        input_dir = "/path/to/data"
        bracken_data = extract_bracken_data(input_dir)
    """
    
    bracken_dir = os.path.join(input_dir, "kraken2_bracken")
    if not os.path.exists(bracken_dir):
        print(f"Error: Bracken directory '{bracken_dir}' does not exist.")
        return {}
    bracken_tsv_files = glob(os.path.join(bracken_dir, "*.bracken.tsv"))
  
    bracken_data = {}
    column_mapping = {
        'name': 'name',
        'taxonomy_id': 'taxid',
        'taxonomy_lvl': 'lvl',
        "added_reads": "ar",
        "new_est_reads": "nr",
        "kraken_assigned_reads":"kar",
        "fraction_total_reads": "frac"}

    for filename in bracken_tsv_files:
        sample_id = os.path.basename(filename).split('.')[0]
        df = pd.read_csv(filename, sep="\t", header=0)
        bracken_data[sample_id] = {}
        # Check if the dataframe is empty
        if df.empty:
            print(f"Warning: The file '{filename}' is empty. Skipping.")
            continue
        df = df.rename(columns=column_mapping)
        # Sort the dataframe by 'fraction_total_reads' in descending order and fetch the first row
        top_row = df.sort_values(by='frac', ascending=False).iloc[0]
        bracken_data[sample_id] = top_row.to_dict()
        bracken_data[sample_id]['frac'] = round(bracken_data[sample_id]['frac'], 2)
    return bracken_data

def extract_plasmidfinder_data(input_dir):
    """
    Extracts and processes PlasmidFinder data from a specified input directory.

    This function searches for TSV files within the "plasmidfinder" subdirectory
    of the given input directory, reads their contents, and processes the data
    into a structured dictionary format. The function also handles cases where
    the directory or files are missing, or the files are empty.

    Args:
        input_dir (str): The path to the input directory containing the 
                         "plasmidfinder" subdirectory.

    Returns:
        dict: A dictionary where keys are sample IDs (derived from the TSV file
              names) and values are lists of dictionaries representing the 
              processed data from each TSV file. If the "plasmidfinder" 
              directory does not exist or no valid data is found, an empty 
              dictionary is returned.

    Notes:
        - The function renames specific columns in the TSV files to match an 
          expected format.
        - Unwanted columns (e.g., 'db') are removed, and missing values are 
          filled with empty strings.
        - If a TSV file is empty, it is skipped with a warning message.

    Example:
        input_dir = "/path/to/input"
        plasmidfinder_data = extract_plasmidfinder_data(input_dir)
        print(plasmidfinder_data)
    """

    plfinder_dir = os.path.join(input_dir, "plasmidfinder")
    if not os.path.exists(plfinder_dir):
        print(f"Error: The PlasmidFinder directory '{plfinder_dir}' does not exist.")
        return {}
    tsv_files = glob(os.path.join(plfinder_dir, "**", "*.tsv"), recursive=True)
    plasmidfinder_data = {}

    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        # Rename columns to match the expected format
        df = df.rename(columns={
            'Database': 'db',
            'Plasmid': 'plasmid',
            'Identity': 'idty',
            'Query / Template length': 'length',
            'Contig': 'contig',
            "Position in contig": "pos",
            "Note": "note",
            "Accession number": "acc"
        })
        # Remove unwanted columns 
        df = df.drop(columns=['db'], errors='ignore')
        df = df.fillna("")
        plasmidfinder_data[sample_id] = df.to_dict(orient='records') 

    return plasmidfinder_data

def get_kleborate_dictionary():
   
    """
    Returns a nested dictionary structure representing the Kleborate report schema.

    The dictionary contains the following main sections:
    - `strain_and_species_identification`: Information about strain and species identification.
        - `strain`: Strain name (default: empty string).
        - `species`: Species name (default: empty string).
        - `species_match`: Species match information (default: empty string).

    - `genome_assembly_metrics`: Metrics related to genome assembly.
        - `contig_count`: Number of contigs (default: 0).
        - `N50`: N50 value of the assembly (default: 0).
        - `largest_contig`: Size of the largest contig (default: 0).
        - `total_size`: Total size of the assembly (default: 0).
        - `ambiguous_bases`: Number of ambiguous bases (default: empty string).
        - `QC_warnings`: Quality control warnings (default: empty string).

    - `sequence_typing`: Sequence typing information.
        - `ST`: Sequence type (default: empty string).
        - `mlst_genes`: Dictionary of MLST genes and their values (default: empty strings).

    - `virulence_factors`: Information about virulence factors.
        - Includes nested dictionaries for `yersiniabactin`, `colibactin`, `aerobactin`, 
          `salmochelin`, and `rmpadc`, each containing subfields for genes, spurious hits, 
          and typing information.
        - `virulence_summary`: Summary of virulence factors, including `virulence_score` 
          and `spurious_virulence_hits`.

    - `antimicrobial_resistance`: Information about antimicrobial resistance.
        - `acquired_resistance_genes`: Dictionary of acquired resistance genes (default: empty strings).
        - `resistance_mutations`: Dictionary of resistance mutations (default: empty strings).
        - `resistance_summary`: Summary of resistance, including scores and counts.

    - `capsule_and_o_antigen_typing`: Capsule and O-antigen typing information.
        - `k_locus`: Information about K locus, including typing and confidence.
        - `o_locus`: Information about O locus, including typing and confidence.

    Returns:
        dict: A dictionary initialized with default values for all Kleborate report fields.
    """
    return {
    "strain_and_species_identification": {
        "strain": "",
        "species": "",
        "species_match": ""
    },
    "genome_assembly_metrics": {
        "contig_count": 0,
        "N50": 0,
        "largest_contig": 0,
        "total_size": 0,
        "ambiguous_bases": "",
        "QC_warnings": ""
    },
    "sequence_typing": {
        "ST": "",
        "mlst_genes": {
            "gapA": "",
            "infB": "",
            "mdh": "",
            "pgi": "",
            "phoE": "",
            "rpoB": "",
            "tonB": ""
        }
    },
    "virulence_factors": {
        "yersiniabactin": {
            "YbST": "",
            "genes": {
                "Yersiniabactin": "",
                "ybtS": "",
                "ybtX": "",
                "ybtQ": "",
                "ybtP": "",
                "ybtA": "",
                "irp2": "",
                "irp1": "",
                "ybtU": "",
                "ybtT": "",
                "ybtE": "",
                "fyuA": ""
            },
            "spurious_ybt_hits": ""
        },
        "colibactin": {
            "CbST": "",
            "genes": {
                "Colibactin": "",
                "clbA": "",
                "clbB": "",
                "clbC": "",
                "clbD": "",
                "clbE": "",
                "clbF": "",
                "clbG": "",
                "clbH": "",
                "clbI": "",
                "clbL": "",
                "clbM": "",
                "clbN": "",
                "clbO": "",
                "clbP": "",
                "clbQ": ""
            },
            "spurious_clb_hits": ""
        },
        "aerobactin": {
            "AbST": "",
            "genes": {
                "Aerobactin": "",
                "iucA": "",
                "iucB": "",
                "iucC": "",
                "iucD": "",
                "iutA": ""
            },
            "spurious_abst_hits": ""
        },
        "salmochelin": {
            "SmST": "",
            "genes": {
                "Salmochelin": "",
                "iroB": "",
                "iroC": "",
                "iroD": "",
                "iroN": ""
            },
            "spurious_smst_hits": ""
        },
        "rmpadc": {
            "RmST": "",
            "genes": {
                "RmpADC": "",
                "rmpA": "",
                "rmpD": "",
                "rmpC": "",
                "rmpA2": ""
            },
            "spurious_rmst_hits": ""
        },
        "virulence_summary": {
            "virulence_score": 0,
            "spurious_virulence_hits": ""
        }
    },
    "antimicrobial_resistance": {
        "acquired_resistance_genes": {
            "AGly_acquired": "",
            "Col_acquired": "",
            "Fcyn_acquired": "",
            "Flq_acquired": "",
            "Gly_acquired": "",
            "MLS_acquired": "",
            "Phe_acquired": "",
            "Rif_acquired": "",
            "Sul_acquired": "",
            "Tet_acquired": "",
            "Tgc_acquired": "",
            "Tmt_acquired": "",
            "Bla_acquired": "",
            "Bla_inhR_acquired": "",
            "Bla_ESBL_acquired": "",
            "Bla_ESBL_inhR_acquired": "",
            "Bla_Carb_acquired": "",
            "Bla_chr": ""
        },
        "resistance_mutations": {
            "SHV_mutations": "",
            "Omp_mutations": "",
            "Col_mutations": "",
            "Flq_mutations": ""
        },
        "resistance_summary": {
            "truncated_resistance_hits": "",
            "spurious_resistance_hits": "",
            "resistance_score": 0,
            "num_resistance_classes": 0,
            "num_resistance_genes": 0
        }
    },
    "capsule_and_o_antigen_typing": {
        "k_locus": {
            "wzi": "",
            "K_locus": "",
            "K_type": "",
            "K_locus_confidence": "",
            "K_locus_problems": "",
            "K_locus_identity": "",
            "K_Missing_expected_genes": ""
        },
        "o_locus": {
            "O_locus": "",
            "O_type": "",
            "O_locus_confidence": "",
            "O_locus_problems": "",
            "O_locus_identity": "",
            "O_Missing_expected_genes": ""
        }
    }
}

def populate_kleborate_dict(row):
    """
    Populate the Kleborate template dictionary with data from a TSV row.
    This function takes a dictionary-like object (`row`) containing data from a TSV file
    and populates a deep copy of the Kleborate template dictionary with the corresponding
    values. If a key is missing in the input row, default values are used.
    Args:
        row (dict): A dictionary-like object containing data from a TSV row. Keys should
                    correspond to the expected fields in the Kleborate template.
    Returns:
        dict: A populated Kleborate dictionary with the data from the input row.
    Notes:
        - The function uses `copy.deepcopy` to ensure the original template dictionary
          remains unmodified.
        - Default values are used for missing keys in the input row:
            - Strings default to an empty string (`""`).
            - Integers default to `0`.
    """
    """Populate the kleborate template dictionary with data from a TSV row."""
    # Create a deep copy of the template to avoid modifying the original
    result = copy.deepcopy(get_kleborate_dictionary())
    
    # Strain and Species Identification
    result["strain_and_species_identification"]["strain"] = row.get("strain", "")
    result["strain_and_species_identification"]["species"] = row.get("species", "")
    result["strain_and_species_identification"]["species_match"] = row.get("species_match", "")
    
    # Genome Assembly Metrics
    result["genome_assembly_metrics"]["contig_count"] = int(row.get("contig_count", 0))
    result["genome_assembly_metrics"]["N50"] = int(row.get("N50", 0))
    result["genome_assembly_metrics"]["largest_contig"] = int(row.get("largest_contig", 0))
    result["genome_assembly_metrics"]["total_size"] = int(row.get("total_size", 0))
    result["genome_assembly_metrics"]["ambiguous_bases"] = row.get("ambiguous_bases", "")
    result["genome_assembly_metrics"]["QC_warnings"] = row.get("QC_warnings", "")
    
    # Sequence Typing
    result["sequence_typing"]["ST"] = row.get("ST", "")
    result["sequence_typing"]["mlst_genes"]["gapA"] = row.get("gapA", "")
    result["sequence_typing"]["mlst_genes"]["infB"] = row.get("infB", "")
    result["sequence_typing"]["mlst_genes"]["mdh"] = row.get("mdh", "")
    result["sequence_typing"]["mlst_genes"]["pgi"] = row.get("pgi", "")
    result["sequence_typing"]["mlst_genes"]["phoE"] = row.get("phoE", "")
    result["sequence_typing"]["mlst_genes"]["rpoB"] = row.get("rpoB", "")
    result["sequence_typing"]["mlst_genes"]["tonB"] = row.get("tonB", "")
    
    # Virulence Factors
    # Yersiniabactin
    result["virulence_factors"]["yersiniabactin"]["YbST"] = row.get("YbST", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["Yersiniabactin"] = row.get("Yersiniabactin", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtS"] = row.get("ybtS", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtX"] = row.get("ybtX", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtQ"] = row.get("ybtQ", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtP"] = row.get("ybtP", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtA"] = row.get("ybtA", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["irp2"] = row.get("irp2", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["irp1"] = row.get("irp1", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtU"] = row.get("ybtU", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtT"] = row.get("ybtT", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["ybtE"] = row.get("ybtE", "")
    result["virulence_factors"]["yersiniabactin"]["genes"]["fyuA"] = row.get("fyuA", "")
    result["virulence_factors"]["yersiniabactin"]["spurious_ybt_hits"] = row.get("spurious_ybt_hits", "")
    
    # Colibactin
    result["virulence_factors"]["colibactin"]["CbST"] = row.get("CbST", "")
    result["virulence_factors"]["colibactin"]["genes"]["Colibactin"] = row.get("Colibactin", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbA"] = row.get("clbA", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbB"] = row.get("clbB", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbC"] = row.get("clbC", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbD"] = row.get("clbD", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbE"] = row.get("clbE", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbF"] = row.get("clbF", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbG"] = row.get("clbG", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbH"] = row.get("clbH", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbI"] = row.get("clbI", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbL"] = row.get("clbL", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbM"] = row.get("clbM", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbN"] = row.get("clbN", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbO"] = row.get("clbO", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbP"] = row.get("clbP", "")
    result["virulence_factors"]["colibactin"]["genes"]["clbQ"] = row.get("clbQ", "")
    result["virulence_factors"]["colibactin"]["spurious_clb_hits"] = row.get("spurious_clb_hits", "")
    
    # Aerobactin
    result["virulence_factors"]["aerobactin"]["AbST"] = row.get("AbST", "")
    result["virulence_factors"]["aerobactin"]["genes"]["Aerobactin"] = row.get("Aerobactin", "")
    result["virulence_factors"]["aerobactin"]["genes"]["iucA"] = row.get("iucA", "")
    result["virulence_factors"]["aerobactin"]["genes"]["iucB"] = row.get("iucB", "")
    result["virulence_factors"]["aerobactin"]["genes"]["iucC"] = row.get("iucC", "")
    result["virulence_factors"]["aerobactin"]["genes"]["iucD"] = row.get("iucD", "")
    result["virulence_factors"]["aerobactin"]["genes"]["iutA"] = row.get("iutA", "")
    result["virulence_factors"]["aerobactin"]["spurious_abst_hits"] = row.get("spurious_abst_hits", "")
    
    # Salmochelin
    result["virulence_factors"]["salmochelin"]["SmST"] = row.get("SmST", "")
    result["virulence_factors"]["salmochelin"]["genes"]["Salmochelin"] = row.get("Salmochelin", "")
    result["virulence_factors"]["salmochelin"]["genes"]["iroB"] = row.get("iroB", "")
    result["virulence_factors"]["salmochelin"]["genes"]["iroC"] = row.get("iroC", "")
    result["virulence_factors"]["salmochelin"]["genes"]["iroD"] = row.get("iroD", "")
    result["virulence_factors"]["salmochelin"]["genes"]["iroN"] = row.get("iroN", "")
    result["virulence_factors"]["salmochelin"]["spurious_smst_hits"] = row.get("spurious_smst_hits", "")
    
    # RmpADC
    result["virulence_factors"]["rmpadc"]["RmST"] = row.get("RmST", "")
    result["virulence_factors"]["rmpadc"]["genes"]["RmpADC"] = row.get("RmpADC", "")
    result["virulence_factors"]["rmpadc"]["genes"]["rmpA"] = row.get("rmpA", "")
    result["virulence_factors"]["rmpadc"]["genes"]["rmpD"] = row.get("rmpD", "")
    result["virulence_factors"]["rmpadc"]["genes"]["rmpC"] = row.get("rmpC", "")
    result["virulence_factors"]["rmpadc"]["genes"]["rmpA2"] = row.get("rmpA2", "")
    result["virulence_factors"]["rmpadc"]["spurious_rmst_hits"] = row.get("spurious_rmst_hits", "")
    
    # Virulence Summary
    result["virulence_factors"]["virulence_summary"]["virulence_score"] = int(row.get("virulence_score", 0))
    result["virulence_factors"]["virulence_summary"]["spurious_virulence_hits"] = row.get("spurious_virulence_hits", "")
    
    # Antimicrobial Resistance
    # Acquired Resistance Genes
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["AGly_acquired"] = row.get("AGly_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Col_acquired"] = row.get("Col_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Fcyn_acquired"] = row.get("Fcyn_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Flq_acquired"] = row.get("Flq_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Gly_acquired"] = row.get("Gly_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["MLS_acquired"] = row.get("MLS_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Phe_acquired"] = row.get("Phe_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Rif_acquired"] = row.get("Rif_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Sul_acquired"] = row.get("Sul_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Tet_acquired"] = row.get("Tet_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Tgc_acquired"] = row.get("Tgc_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Tmt_acquired"] = row.get("Tmt_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Bla_acquired"] = row.get("Bla_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Bla_inhR_acquired"] = row.get("Bla_inhR_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Bla_ESBL_acquired"] = row.get("Bla_ESBL_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Bla_ESBL_inhR_acquired"] = row.get("Bla_ESBL_inhR_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Bla_Carb_acquired"] = row.get("Bla_Carb_acquired", "")
    result["antimicrobial_resistance"]["acquired_resistance_genes"]["Bla_chr"] = row.get("Bla_chr", "")
    
    # Resistance Mutations
    result["antimicrobial_resistance"]["resistance_mutations"]["SHV_mutations"] = row.get("SHV_mutations", "")
    result["antimicrobial_resistance"]["resistance_mutations"]["Omp_mutations"] = row.get("Omp_mutations", "")
    result["antimicrobial_resistance"]["resistance_mutations"]["Col_mutations"] = row.get("Col_mutations", "")
    result["antimicrobial_resistance"]["resistance_mutations"]["Flq_mutations"] = row.get("Flq_mutations", "")
    
    # Resistance Summary
    result["antimicrobial_resistance"]["resistance_summary"]["truncated_resistance_hits"] = row.get("truncated_resistance_hits", "")
    result["antimicrobial_resistance"]["resistance_summary"]["spurious_resistance_hits"] = row.get("spurious_resistance_hits", "")
    result["antimicrobial_resistance"]["resistance_summary"]["resistance_score"] = int(row.get("resistance_score", 0))
    result["antimicrobial_resistance"]["resistance_summary"]["num_resistance_classes"] = int(row.get("num_resistance_classes", 0))
    result["antimicrobial_resistance"]["resistance_summary"]["num_resistance_genes"] = int(row.get("num_resistance_genes", 0))
    
    # Capsule and O-Antigen Typing
    # K-Locus
    result["capsule_and_o_antigen_typing"]["k_locus"]["wzi"] = row.get("wzi", "")
    result["capsule_and_o_antigen_typing"]["k_locus"]["K_locus"] = row.get("K_locus", "")
    result["capsule_and_o_antigen_typing"]["k_locus"]["K_type"] = row.get("K_type", "")
    result["capsule_and_o_antigen_typing"]["k_locus"]["K_locus_confidence"] = row.get("K_locus_confidence", "")
    result["capsule_and_o_antigen_typing"]["k_locus"]["K_locus_problems"] = row.get("K_locus_problems", "")
    result["capsule_and_o_antigen_typing"]["k_locus"]["K_locus_identity"] = row.get("K_locus_identity", "")
    result["capsule_and_o_antigen_typing"]["k_locus"]["K_Missing_expected_genes"] = row.get("K_Missing_expected_genes", "")
    
    # O-Locus
    result["capsule_and_o_antigen_typing"]["o_locus"]["O_locus"] = row.get("O_locus", "")
    result["capsule_and_o_antigen_typing"]["o_locus"]["O_type"] = row.get("O_type", "")
    result["capsule_and_o_antigen_typing"]["o_locus"]["O_locus_confidence"] = row.get("O_locus_confidence", "")
    result["capsule_and_o_antigen_typing"]["o_locus"]["O_locus_problems"] = row.get("O_locus_problems", "")
    result["capsule_and_o_antigen_typing"]["o_locus"]["O_locus_identity"] = row.get("O_locus_identity", "")
    result["capsule_and_o_antigen_typing"]["o_locus"]["O_Missing_expected_genes"] = row.get("O_Missing_expected_genes", "")
    
    return result

def extract_kleborate_data(input_dir):
    """
    Extracts Kleborate data from TSV files located in the 'kleborate' subdirectory of the given input directory.
    This function searches for TSV files in the 'kleborate' directory, reads their contents, and processes
    the data into a dictionary where each key is a sample ID (derived from the filename) and the value is
    the processed data for that sample.
    Args:
        input_dir (str): The path to the input directory containing the 'kleborate' subdirectory.
    Returns:
        dict: A dictionary where keys are sample IDs and values are the processed Kleborate data.
              Returns an empty dictionary if the 'kleborate' directory does not exist or contains no TSV files.
    Raises:
        FileNotFoundError: If the specified 'kleborate' directory does not exist.
        csv.Error: If there is an error reading a TSV file.
    Notes:
        - The function assumes that TSV files use tab ('\t') as the delimiter.
        - The `populate_kleborate_dict` function is used to process each row of the TSV file.
    """


    kleborate_dir = os.path.join(input_dir, "kleborate")
    if not os.path.exists(kleborate_dir):
        print(f"Error: The Kleborate directory '{kleborate_dir}' does not exist.")
        return {}
    
    tsv_files = glob(os.path.join(kleborate_dir, "*.tsv"))
    if not tsv_files:
        print(f"Error: No TSV files found in '{kleborate_dir}'.")
        return {}
    kleborate_data = {}
    for tsv_file in tsv_files: 
        sample_id = os.path.basename(tsv_file).split('.')[0]
        with open(tsv_file, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                kleborate_data[sample_id] = populate_kleborate_dict(row)    
    return kleborate_data

def extract_vfdb_data(input_dir):
    
    """
    This function reads TSV files from the `diamond/VFDB` subdirectory within the given `input_dir`.
    It processes each file to extract virulence factor information, filters the data based on 
    specific criteria, and returns the results in a structured dictionary format.
    The filtering criteria include:
    - Percentage identity (`pident`) greater than 70.
    - Subject coverage (`scov`) greater than 50%.
    - Query coverage (`qcov`) greater than 50%.
    Additionally, the function extracts gene names and notes from the `stitle` column of the TSV files.
        input_dir (str): Path to the input directory containing the `diamond/VFDB` folder.
        dict: A dictionary where:
            - Keys are sample IDs (derived from file names).
            - Values are lists of dictionaries representing the processed data from each TSV file.
    Notes:
        - If the `diamond/VFDB` directory does not exist, the function prints an error message and returns an empty dictionary.
        - If a TSV file is empty, the function skips processing it and prints a warning message.
        - The function drops unnecessary columns from the processed data, including `stitle`, `evalue`, `qstart`, `qend`, 
          `sstart`, `send`, `slen`, and `qlen`.
    Raises:
        None
    """
    
    vfdb_dir = os.path.join(input_dir, "diamond/VFDB")
    if not os.path.exists(vfdb_dir):
        print(f"Error: VFDB directory '{vfdb_dir}' does not exist.")
        return {}
    
    tsv_files = glob(os.path.join(vfdb_dir, "*.tsv"))
    vfdb_data = {}

    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        df["scov"] = 100. * round(df["length"] / df["slen"], 2)
        df["qcov"] = 100. * round(3 * df["length"] / df["qlen"], 2)
        df = df[(df["pident"] > 70) & (df["scov"] > 50) & (df["qcov"] > 50)]
        df["gene_name"] = df["stitle"].str.extract(r'\)\s*\((.*?)\)', expand=False)
        df["note"] = df["stitle"].str.replace(r'^.*?\)\s*\((.*?)\)\s*', '', regex=True)
        df = df.drop(columns=["stitle", "evalue", "qstart", "qend", "sstart", "send","slen","qlen"], errors='ignore')
        vfdb_data[sample_id] = df.to_dict(orient='records')

    return vfdb_data

def extract_rgi_data(input_dir):

    rgi_dir = os.path.join(input_dir, "rgi")
    if not os.path.exists(rgi_dir):
        print(f"Error: The RGI directory '{rgi_dir}' does not exist.")
        return {}

    tsv_files = glob(os.path.join(rgi_dir, "**", "*.rgi.txt"), recursive=True)
    if not tsv_files:
        print(f"Error: No RGI TSV files found in '{rgi_dir}'.")
        return {}

    rgi_data = {}

    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]

        try:
            df = pd.read_csv(tsv_file, sep="\t", header=0)
        except Exception as e:
            print(f"Error reading file '{tsv_file}': {e}")
            continue

        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue

        # Drop unnecessary columns
        columns_to_drop = [
            "Contig", "Start", "Stop", "Orientation", "Predicted_Protein",
            "CARD_Protein_Sequence ID", "Model_ID", "Nudged", "Note",
            "Hit_Start", "Hit_End", "Antibiotic","CARD_Protein_Sequence","Predicted_DNA"
        ]
        df = df.drop(columns=columns_to_drop, errors='ignore')

        # Filter out rows based on Best_Identities greater than 70
        df = df[df["Best_Identities"] >= 70]

        # Rename columns for consistency
        column_renames = {
            "Drug Class": "Drug_Class",
            "Resistance Mechanism": "Resistance_Mechanism",
            "AMR Gene Family": "AMR_Gene_Family",
            "Percentage Length of Reference Sequence": "scov"
        }
        df = df.rename(columns=column_renames)

        # Clean up specific columns
        df["ORF_ID"] = df["ORF_ID"].str.split(" ").str[0]
        df["SNPs_in_Best_Hit_ARO"] = df["SNPs_in_Best_Hit_ARO"].replace("n/a", "")
        df["Other_SNPs"] = df["Other_SNPs"].replace("n/a", "")

        # Convert the dataframe to a dictionary
        rgi_data[sample_id] = df.to_dict(orient='records')

    return rgi_data

def extract_amrfinder_data(input_dir):
    """
    Extracts and processes AMRFinder data from a specified input directory.
    This function reads TSV files generated by AMRFinder from a subdirectory named "amrfinder"
    within the given input directory. It processes the data to extract information about 
    antimicrobial resistance (AMR), virulence factors (VF), and stress response elements, 
    organizing the results into a structured dictionary.
    Args:
        input_dir (str): The path to the input directory containing the "amrfinder" subdirectory.
    Returns:
        dict: A dictionary where each key is a sample ID (derived from the TSV filenames) and 
              the value is another dictionary with keys 'amr', 'vf', and 'stress', each containing 
              a list of records (as dictionaries) for the respective element type. If the 
              "amrfinder" directory or TSV files are missing, or if the files are empty, 
              an empty dictionary is returned.
    Notes:
        - The function expects the TSV files to have specific column names as generated by AMRFinder.
        - Columns are renamed for consistency, and unnecessary columns are dropped.
        - Warnings are printed if files are empty or if no TSV files are found.
    Example:
        input_dir/
        ├── amrfinder/
        │   ├── sample1.amrfinder.tsv
        │   ├── sample2.amrfinder.tsv
        amrfinder_data = extract_amrfinder_data(input_dir)
        print(amrfinder_data)
    """
    amrfinder_dir = os.path.join(input_dir, "amrfinder")
    if not os.path.exists(amrfinder_dir):
        print(f"Error: The AMRFinder directory '{amrfinder_dir}' does not exist.")
        return {}
    
    tsv_files = glob(os.path.join(amrfinder_dir, "*.amrfinder.tsv"))      
    if not tsv_files:  
        print(f"Error: No AMRFinder TSV files found in '{amrfinder_dir}'.")
        return {}
    amrfinder_data = {}
    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        amrfinder_data[sample_id] = {'amr': [], 'vf': [], 'stress': []}
        # Drop unnecessary columns
        df = df.drop(columns=["Protein identifier", "HMM id", "HMM description",], errors='ignore')
       
        # Rename columns for consistency
        column_renames = {
            "Contig id": "cid",
            "Start": "sidx",
            "Stop": "eidx",
            "Strand": "snd",
            "Cut_Off": "cutoff",
            "Gene symbol": "gene_name",
            "Sequence name": "product",
            "Scope": "scope",
            "Element type": "element_type",
            "Element subtype": "element_subtype",
            "Class": "class",
            "Subclass": "subclass",
            "Method": "method",
            "Target length": "qlen",
            "Reference sequence length": "slen",
            "% Coverage of reference sequence": "scov",
            "% Identity to reference sequence": "idty",
            "Alignment length": "length",
            "Accession of closest sequence": "acc",
            "Name of closest sequence": "dscr"
        }
        df = df.rename(columns=column_renames)
        
        df_amr = df[df["element_type"] == "AMR"]
        df_vf = df[df["element_type"] == "VIRULENCE"]
        df_stress = df[df["element_type"] == "STRESS"]
        amrfinder_data[sample_id]['amr'] = df_amr.to_dict(orient='records')
        amrfinder_data[sample_id]['vf'] = df_vf.to_dict(orient='records')
        amrfinder_data[sample_id]['stress'] = df_stress.to_dict(orient='records')

    return amrfinder_data

def extract_amrfinder_mutations_data(input_dir):
    """
   
    """
    amrfinder_dir = os.path.join(input_dir, "amrfinder")
    if not os.path.exists(amrfinder_dir):
        print(f"Error: The AMRFinder directory '{amrfinder_dir}' does not exist.")
        return {}
    
    tsv_files = glob(os.path.join(amrfinder_dir, "*.mutations.tsv"))      
    if not tsv_files:  
        print(f"Error: No AMRFinder TSV files found in '{amrfinder_dir}'.")
        return {}
    amrfinder_data = {}
    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        amrfinder_data[sample_id] = {'amr': [], 'vf': [], 'stress': []}
        # Drop unnecessary columns
        df = df.drop(columns=["Protein identifier", "HMM id", "HMM description",], errors='ignore')
       
        # Rename columns for consistency
        column_renames = {
            "Contig id": "cid",
            "Start": "sidx",
            "Stop": "eidx",
            "Strand": "snd",
            "Cut_Off": "cutoff",
            "Gene symbol": "gene_name",
            "Sequence name": "product",
            "Scope": "scope",
            "Element type": "element_type",
            "Element subtype": "element_subtype",
            "Class": "class",
            "Subclass": "subclass",
            "Method": "method",
            "Target length": "qlen",
            "Reference sequence length": "slen",
            "% Coverage of reference sequence": "scov",
            "% Identity to reference sequence": "idty",
            "Alignment length": "length",
            "Accession of closest sequence": "acc",
            "Name of closest sequence": "dscr"
        }
        df = df.rename(columns=column_renames)
        
       
        amrfinder_data[sample_id] = df.to_dict(orient='records')

    return amrfinder_data

def extract_snippy_data(input_dir):
    """
    Extracts and processes mutation data from Snippy CSV files within a specified input directory.
    This function searches for CSV files in the "snippy" subdirectory of the given input directory,
    filters out synonymous variants and non-CDS features, and organizes the mutation data into a
    dictionary format.
    Args:
        input_dir (str): The path to the input directory containing the "snippy" subdirectory.
    Returns:
        dict: A dictionary where keys are sample IDs (derived from CSV filenames) and values are
              dictionaries containing:
              - "reference" (str): The reference genome name (from the "CHROM" column).
              - "mutations" (list): A list of dictionaries representing mutations, with columns
                from the CSV file as keys (excluding dropped columns).
    Notes:
        - If the "snippy" directory does not exist or no CSV files are found, an empty dictionary
          is returned.
        - If a CSV file cannot be read or is empty, it is skipped with a warning.
        - The following columns are dropped from the mutation data if present: "CHROM", "TYPE",
          "FTYPE", "STRAND", "NT_POS", "AA_POS", "EVIDENCE".
        - Only rows where "EFFECT" does not contain "synonymous_variant" and "FTYPE" is "CDS" are
          included in the mutation data.
    Raises:
        None: Any exceptions encountered while reading CSV files are caught and logged as errors.
    """

    snippy_dir = os.path.join(input_dir, "snippy")
    if not os.path.exists(snippy_dir):
        print(f"Error: The snippy directory '{snippy_dir}' does not exist.")
        return {}
    
    csv_files = glob(os.path.join(snippy_dir, "**", "*.csv"), recursive=True)
    if not csv_files:
        print(f"Error: No snippy CSV files found in '{snippy_dir}'.")
        return {}
    snippy_data = {}
    for csv_file in csv_files:
        sample_id = os.path.basename(csv_file).split('.')[0]
        try:
            df = pd.read_csv(csv_file, header=0)
        except Exception as e:
            print(f"Error reading file '{csv_file}': {e}")
            continue

        if df.empty:
            print(f"Warning: The file '{csv_file}' is empty. Skipping.")
            continue
        
        snippy_data[sample_id]= {"reference": "", "mutations": []}
        # Assign the value of the CHROM column to the reference field
        if "CHROM" in df.columns:
            snippy_data[sample_id]["reference"] = df["CHROM"].iloc[0]

        
        
        # Filter rows where EFFECT does contain "synonymous_variant" and FTYPE is 'CDS'
        if "EFFECT" in df.columns and "FTYPE" in df.columns:
            df = df[~df["EFFECT"].str.contains("synonymous_variant", na=False)]
        
        # Drop unnecessary columns
        columns_to_drop = ["CHROM","TYPE","FTYPE","STRAND","NT_POS","AA_POS","EVIDENCE"]
        df = df.drop(columns=columns_to_drop, errors='ignore')

        snippy_data[sample_id]["mutations"] = df.to_dict(orient='records')

    return snippy_data

def get_quast_dictionary():
    return {
    "contig_count_metrics": {
        "# contigs (>= 0 bp)": 0,
        "# contigs (>= 1000 bp)": 0,
        "# contigs (>= 5000 bp)": 0,
        "# contigs (>= 10000 bp)": 0,
        "# contigs (>= 25000 bp)": 0,
        "# contigs (>= 50000 bp)": 0,
        "# contigs": 0
    },
    "contig_length_metrics": {
        "Total length (>= 0 bp)": 0,
        "Total length (>= 1000 bp)": 0,
        "Total length (>= 5000 bp)": 0,
        "Total length (>= 10000 bp)": 0,
        "Total length (>= 25000 bp)": 0,
        "Total length (>= 50000 bp)": 0,
        "Total length": 0,
        "Largest contig": 0
    },
    "assembly_statistics": {
        "N50": 0,
        "NG50": 0,
        "N90": 0,
        "NG90": 0,
        "auN": 0.0,
        "auNG": 0.0,
        "L50": 0,
        "LG50": 0,
        "L90": 0,
        "LG90": 0
    },
    "reference_based_metrics": {
        "Reference length": 0,
        "GC (%)": 0.0,
        "Reference GC (%)": 0.0,
        "Genome fraction (%)": 0.0,
        "Duplication ratio": 0.0
    },
    "misassembly_metrics": {
        "# misassemblies": 0,
        "# misassembled contigs": 0,
        "Misassembled contigs length": 0,
        "# local misassemblies": 0,
        "# scaffold gap ext. mis.": 0,
        "# scaffold gap loc. mis.": 0,
        "# unaligned mis. contigs": 0,
        "# unaligned contigs": 0,
        "Unaligned length": 0
    },
    "alignment_based_statistics": {
        "Largest alignment": 0,
        "Total aligned length": 0,
        "NA50": 0,
        "NGA50": 0,
        "NA90": 0,
        "NGA90": 0,
        "auNA": 0.0,
        "auNGA": 0.0,
        "LA50": 0,
        "LGA50": 0,
        "LA90": 0,
        "LGA90": 0
    },
    "error_metrics": {
        "# N's per 100 kbp": 0.0,
        "# mismatches per 100 kbp": 0.0,
        "# indels per 100 kbp": 0.0
    }
}

def populate_quast_dict(row):

    result = copy.deepcopy(get_quast_dictionary())
    
    # Contig Count Metrics
    result["contig_count_metrics"]["# contigs (>= 0 bp)"] = row.get("# contigs (>= 0 bp)", "")
    result["contig_count_metrics"]["# contigs (>= 1000 bp)"] = row.get("# contigs (>= 1000 bp)", "")
    result["contig_count_metrics"]["# contigs (>= 5000 bp)"] = row.get("# contigs (>= 5000 bp)", "")
    result["contig_count_metrics"]["# contigs (>= 10000 bp)"] = row.get("# contigs (>= 10000 bp)", "")
    result["contig_count_metrics"]["# contigs (>= 25000 bp)"] = row.get("# contigs (>= 25000 bp)", "")
    result["contig_count_metrics"]["# contigs (>= 50000 bp)"] = row.get("# contigs (>= 50000 bp)", "")
    result["contig_count_metrics"]["# contigs"] = row.get("# contigs", "")
    
    # Contig Length Metrics
    result["contig_length_metrics"]["Total length (>= 0 bp)"] = row.get("Total length (>= 0 bp)", "")
    result["contig_length_metrics"]["Total length (>= 1000 bp)"] = row.get("Total length (>= 1000 bp)", "")
    result["contig_length_metrics"]["Total length (>= 5000 bp)"] = row.get("Total length (>= 5000 bp)", "")
    result["contig_length_metrics"]["Total length (>= 10000 bp)"] = row.get("Total length (>= 10000 bp)", "")
    result["contig_length_metrics"]["Total length (>= 25000 bp)"] = row.get("Total length (>= 25000 bp)", "")
    result["contig_length_metrics"]["Total length (>= 50000 bp)"] = row.get("Total length (>= 50000 bp)", "")
    result["contig_length_metrics"]["Total length"] = row.get("Total length", "")
    result["contig_length_metrics"]["Largest contig"] = row.get("Largest contig", "")
    
    # Assembly Statistics
    result["assembly_statistics"]["N50"] = row.get("N50", "")
    result["assembly_statistics"]["NG50"] = row.get("NG50", "")
    result["assembly_statistics"]["N90"] = row.get("N90", "")
    result["assembly_statistics"]["NG90"] = row.get("NG90", "")
    result["assembly_statistics"]["auN"] = row.get("auN", "")
    result["assembly_statistics"]["auNG"] = row.get("auNG", "")
    result["assembly_statistics"]["L50"] = row.get("L50", "")
    result["assembly_statistics"]["LG50"] = row.get("LG50", "")
    result["assembly_statistics"]["L90"] = row.get("L90", "")
    result["assembly_statistics"]["LG90"] = row.get("LG90", "")
    
    # Reference-Based Metrics
    result["reference_based_metrics"]["Reference length"] = row.get("Reference length", "")
    result["reference_based_metrics"]["GC (%)"] = row.get("GC (%)", "")
    result["reference_based_metrics"]["Reference GC (%)"] = row.get("Reference GC (%)", "")
    result["reference_based_metrics"]["Genome fraction (%)"] = row.get("Genome fraction (%)", "")
    result["reference_based_metrics"]["Duplication ratio"] = row.get("Duplication ratio", "")
    
    # Misassembly Metrics
    result["misassembly_metrics"]["# misassemblies"] = row.get("# misassemblies", "")
    result["misassembly_metrics"]["# misassembled contigs"] = row.get("# misassembled contigs", "")
    result["misassembly_metrics"]["Misassembled contigs length"] = row.get("Misassembled contigs length", "")
    result["misassembly_metrics"]["# local misassemblies"] = row.get("# local misassemblies", "")
    result["misassembly_metrics"]["# scaffold gap ext. mis."] = row.get("# scaffold gap ext. mis.", "")
    result["misassembly_metrics"]["# scaffold gap loc. mis."] = row.get("# scaffold gap loc. mis.", "")
    result["misassembly_metrics"]["# unaligned mis. contigs"] = row.get("# unaligned mis. contigs", "")
    result["misassembly_metrics"]["# unaligned contigs"] = row.get("# unaligned contigs", "")
    result["misassembly_metrics"]["Unaligned length"] = row.get("Unaligned length", "")
    
    # Alignment-Based Statistics
    result["alignment_based_statistics"]["Largest alignment"] = row.get("Largest alignment", "")
    result["alignment_based_statistics"]["Total aligned length"] = row.get("Total aligned length", "")
    result["alignment_based_statistics"]["NA50"] = row.get("NA50", "")
    result["alignment_based_statistics"]["NGA50"] = row.get("NGA50", "")
    result["alignment_based_statistics"]["NA90"] = row.get("NA90", "")
    result["alignment_based_statistics"]["NGA90"] = row.get("NGA90", "")
    result["alignment_based_statistics"]["auNA"] = row.get("auNA", "")
    result["alignment_based_statistics"]["auNGA"] = row.get("auNGA", "")
    result["alignment_based_statistics"]["LA50"] = row.get("LA50", "")
    result["alignment_based_statistics"]["LGA50"] = row.get("LGA50", "")
    result["alignment_based_statistics"]["LA90"] = row.get("LA90", "")
    result["alignment_based_statistics"]["LGA90"] = row.get("LGA90", "")
    
    # Error Metrics
    result["error_metrics"]["# N's per 100 kbp"] = row.get("# N's per 100 kbp", "")
    result["error_metrics"]["# mismatches per 100 kbp"] = row.get("# mismatches per 100 kbp", "")
    result["error_metrics"]["# indels per 100 kbp"] = row.get("# indels per 100 kbp", "")
    
    return result

def extract_quast_data(input_dir):

    """
    Extracts QUAST data from TSV files located in the 'quast' subdirectory of the given input directory.
    This function searches for TSV files in the 'quast' directory, reads their contents, and processes
    the data into a dictionary where each key is a sample ID (derived from the filename) and the value is
    the processed data for that sample.
    Args:
        input_dir (str): The path to the input directory containing the 'quast' subdirectory.
    Returns:
        dict: A dictionary where keys are sample IDs and values are the processed QUAST data.
              Returns an empty dictionary if the 'quast' directory does not exist or contains no TSV files.
    Raises:
        FileNotFoundError: If the specified 'quast' directory does not exist.
        csv.Error: If there is an error reading a TSV file.
    Notes:
        - The function assumes that TSV files use tab ('\t') as the delimiter.
        - The `populate_quast_dict` function is used to process each row of the TSV file.
    """

    quast_dir = os.path.join(input_dir, "quast")
    if not os.path.exists(quast_dir):
        print(f"Error: The QUAST directory '{quast_dir}' does not exist.")
        return {}
    
    tsv_files = glob(os.path.join(quast_dir, "**", "transposed_report.tsv"), recursive=True)
    if not tsv_files:
        print(f"Error: No QUAST TSV files found in '{quast_dir}'.")
        return {}
    quast_data = {}
    for tsv_file in tsv_files:
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue

        if "Assembly" in df.columns:
            sample_id = str(df["Assembly"].iloc[0]).split('.')[0]
        else:
            sample_id = os.path.basename(os.path.dirname(tsv_file))

        for _, row in df.iterrows():
            # Convert row to dictionary, ensuring missing values are handled
            row_dict = row.to_dict()
            quast_data[sample_id] = populate_quast_dict(row_dict)

    return quast_data

def extract_fatqc_summary(tsv_file):
    """Read a TSV file with FastQC results and convert to a nested dictionary."""
    # Read TSV into a pandas DataFrame
    df = pd.read_csv(tsv_file, sep='\t', header=None, names=['Status', 'Metric', 'Sample'])
    
    # Initialize the result dictionary
    result = {}
    
    # Process each row
    # Extract the sample name from the first row
    sample_full = df.iloc[0]['Sample']
    if sample_full.endswith('_1.gz') or sample_full.endswith('_R1.gz'):
        sample_name = '_'.join(sample_full.split("_")[:-1])
        read_pair = 'fr'
    elif sample_full.endswith('_2.gz') or sample_full.endswith('_R2.gz'):
        sample_name = '_'.join(sample_full.split("_")[:-1])
        read_pair = 'rv'
    else:
        print("Error: Sample name format not recognized.")
        return result  # Return empty result if format is not recognized

    # Initialize sample dictionary
    result = {
        "sample_name": sample_name,
        "read_pair": read_pair,
        "metrics": {}
    }

    for _, row in df.iterrows():
        status = row['Status']
        metric = row['Metric']
        # Add metric and status to the metrics dictionary
        result["metrics"][metric] = status

    return result

def extract_fastqc_data(input_dir, img_dir):
    # Remove the temporary directory
    def clean_tmp_dir(directory):
        """Remove the content of the specified temporary directory."""
        for root, dirs, files in os.walk(directory, topdown=False):
            for file in files:
                os.remove(os.path.join(root, file))
            for dir in dirs:
                os.rmdir(os.path.join(root, dir))


    # Find all .zip files in the input_dir/fastqc directory
    fastqc_dir = os.path.join(input_dir, "fastqc")
    zip_files = glob(os.path.join(fastqc_dir ,"*fastqc.zip"))

    if not zip_files:
        print(f"Error: No FastQC zip files found in '{fastqc_dir}'.")
        return {}

        # Create img_dir if it does not exist
    os.makedirs(img_dir, exist_ok=True)

    # Create a temporary directory for extraction
    tmp_dir = os.path.join(img_dir, "tmp_fastqc")
    os.makedirs(tmp_dir, exist_ok=True)
   
    # Clean the temporary directory if it exists
    clean_tmp_dir(tmp_dir)

    fastqc_data = {}

    for zip_file in zip_files:
        read_name = os.path.basename(zip_file).split('.')[0]

        # Extract the zip file to the temporary directory
        with zipfile.ZipFile(zip_file, 'r') as zip_ref:
            zip_ref.extractall(tmp_dir)

        # Locate the summary.txt file
        summary_files = glob(os.path.join(tmp_dir, '**', "summary.txt"), recursive=True)

        if len(summary_files) == 0:
            print(f"Warning: summary.txt not found for '{read_name}'. Skipping.")
            continue

        # Process the summary.txt file
        fastqc_read = extract_fatqc_summary(summary_files[0])
        if fastqc_read["sample_name"] not in fastqc_data:
            fastqc_data[fastqc_read["sample_name"]] = {}
        if fastqc_read["read_pair"] not in fastqc_data[fastqc_read["sample_name"]]:
            fastqc_data[fastqc_read["sample_name"]][fastqc_read["read_pair"]] = {}

        fastqc_data[fastqc_read["sample_name"]][fastqc_read["read_pair"]] = fastqc_read["metrics"]

        # Locate the images and copy them to img_dir
        images = glob(os.path.join(tmp_dir, "**","**", "per_base_quality.svg"))
        if len(images) > 0:  
            image_name = f"{fastqc_read['sample_name']}_{fastqc_read['read_pair']}.svg"
            dest_image_path = os.path.join(img_dir, image_name)
            shutil.copy(images[0], dest_image_path)
        else:
            print(f"Warning: No images found for '{read_name}'. Skipping.")

        clean_tmp_dir(tmp_dir)

    # Remove the temporary directory if it exists
    if os.path.exists(tmp_dir):
            os.rmdir(tmp_dir)


    return fastqc_data

def generate_report(input_dir, data_dir, img_dir):
    """
    Generate a report from the collected data.
    This is a placeholder function and should be implemented based on your requirements.
    """
    data = {}
    print("=" * 50)
    print("🔍 Extracting MLST data...")
    data["mlst"] = extract_mlst_data(input_dir)
    print("✅ MLST data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting Bracken data...")
    data["bracken"] = extract_bracken_data(input_dir)
    print("✅ Bracken data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting PlasmidFinder data...")
    data["plasmidfinder"] = extract_plasmidfinder_data(input_dir)
    print("✅ PlasmidFinder data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting Kleborate data...")
    data["kleborate"] = extract_kleborate_data(input_dir)
    print("✅ Kleborate data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting VFDB data...")
    data["vfdb"] = extract_vfdb_data(input_dir)
    print("✅ VFDB data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting RGI data...")
    data["rgi"] = extract_rgi_data(input_dir)
    print("✅ RGI data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting AMRFinder data...")
    data["amrfinder"] = extract_amrfinder_data(input_dir)
    print("✅ AMRFinder data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting AMRFinder mutations data...")
    data["amrfinder_mutations"] = extract_amrfinder_mutations_data(input_dir)
    print("✅ AMRFinder mutations data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting QUAST data...")
    data["quast"] = extract_quast_data(input_dir)
    print("✅ QUAST data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting FastQC data...")
    data["fastqc"] = extract_fastqc_data(input_dir, os.path.join(img_dir, "qc"))
    print("✅ FastQC data extraction completed.")
    print("=" * 50)

    print("🔍 Extracting Snippy data...")
    data["snippy"] = extract_snippy_data(input_dir)
    print("✅ Snippy data extraction completed.")
    print("=" * 50)

    # Save the extracted data in JSON format
    output_file = os.path.join(data_dir, "data.json")
    with open(output_file, 'w') as json_file:
        json.dump(data, json_file)
    print(f"Report saved to '{output_file}'.")

    # Convert the data dictionary to a JSON string
    json_string = json.dumps(data)

    # Save the JSON string to a js file
    js_string_file = os.path.join(data_dir, "data.js")
    with open(js_string_file, 'w') as js_file:
        js_file.write(f"const data = {json_string};")
    print(f"Data as JSON string saved to '{js_string_file}'.")

def main(input_dir, output_dir):
    if not os.path.exists(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory '{output_dir}'.")

    # Add your processing logic here
    print(f"Processing files from '{input_dir}' and saving results to '{output_dir}'.")
    # Create subdirectories for images and data
    img_dir = os.path.join(output_dir, "img")
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(img_dir, exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)

    # Generate the report
    generate_report(input_dir, data_dir, img_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process input and output directories.")
    parser.add_argument("--input_dir", "-i", type=str, required=True, help="Path to the input directory.")
    parser.add_argument("--output_dir", "-o", type=str, required=True, help="Path to the output directory.")
    args = parser.parse_args()

    main(args.input_dir, args.output_dir)