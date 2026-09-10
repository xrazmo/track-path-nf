
import os
import sqlite3 as sql
import argparse
import pandas as pd
from glob import glob
from datetime import datetime
import hashlib
import json
from Bio import SeqIO
from vcf_parser import process_vcf_files_parallel

ID_MAPPER = {}

def create_table_if_not_exists(conn, table_name, df):
    """Create a table based on dataframe structure if it doesn't exist already."""
    # Get column names and types from dataframe
    columns = []
    
    for col_name in df.columns:
        # Map pandas dtypes to SQLite types
        if pd.api.types.is_integer_dtype(df[col_name]):
            col_type = "INTEGER"
        elif pd.api.types.is_float_dtype(df[col_name]):
            col_type = "REAL"
        elif pd.api.types.is_datetime64_any_dtype(df[col_name]):
            col_type = "TEXT"  # Store datetime as text in ISO format
        else:
            col_type = "TEXT"
        
        columns.append(f'"{col_name}" {col_type}')
    
    # Add timestamp and record_hash columns
    columns.append('"import_timestamp" TEXT')
    columns.append('"record_hash" TEXT')
    
    # Create table with columns
    create_table_sql = f'''
    CREATE TABLE IF NOT EXISTS "{table_name}" (
        {", ".join(columns)},
        PRIMARY KEY ("record_hash")
    )
    '''
    
    conn.execute(create_table_sql)
    conn.commit()
    print(f"Table '{table_name}' structure created or verified.")

def calculate_record_hash(row,include_index=True):
    """
    Calculate a more robust hash for a record with extremely low collision probability.
    
    This function creates a more robust hash by:
    1. Using SHA256 instead of MD5
    2. Including data types in the hash calculation
    3. Optional inclusion of the row index for absolute uniqueness
    
    Args:
        row (pandas.Series): Row of data to hash
        include_index (bool): Whether to include row index in hash calculation
    
    Returns:
        str: Hex digest of the hash
    """
    # Create a list to store the values with their types
    typed_values = []
    
    # Process each value in the row
    for column, value in row.items():
        if column == 'record_hash' or column == 'exists' or column == 'import_timestamp':
            continue  # Skip metadata columns
        
        if column == 'nt_seq' or column == 'aa_seq':
            continue  # Skip sequence columns
        
        # Add type information to each value
        value_type = type(value).__name__
        
        if pd.isna(value):
            typed_values.append(f"NULL:{column}")
        elif isinstance(value, (int, float, bool)):
            typed_values.append(f"{value_type}:{column}:{value}")
        else:
            # For strings and other types, convert to string
            typed_values.append(f"{value_type}:{column}:{str(value)}")
    
    # Add index if requested
    if include_index and hasattr(row, 'name'):
        typed_values.append(f"index:{row.name}")
    
    # Sort to ensure consistent ordering
    typed_values.sort()
    
    # Create a single string and hash it with SHA-256
    row_str = json.dumps(typed_values)
    return hashlib.sha256(row_str.encode()).hexdigest()

def insert_data(conn, table_name, df, add_timestamp=True):
    """
    Insert dataframe rows into SQLite with deduplication using
    calculate_record_hash().
    """

    if df.empty:
        return 0

    df = df.copy()

    # ------------------------------------------------------------
    # Add import timestamp (once per batch)
    # ------------------------------------------------------------
    if add_timestamp:
        df["import_timestamp"] = datetime.now().isoformat()

    # ------------------------------------------------------------
    # Compute record_hash using custom function
    # ------------------------------------------------------------
    if "record_hash" not in df.columns:
        df["record_hash"] = df.apply(
            calculate_record_hash,
            axis=1,
            include_index=False  # CRITICAL: avoid pandas index instability
        )

    # ------------------------------------------------------------
    # Find which hashes are new (SQLite-side)
    # ------------------------------------------------------------
    cursor = conn.cursor()
    temp_table = f"temp_hashes_{table_name}"

    cursor.execute(f"""
        CREATE TEMPORARY TABLE {temp_table} (
            record_hash TEXT PRIMARY KEY
        )
    """)

    cursor.executemany(
        f"INSERT OR IGNORE INTO {temp_table} (record_hash) VALUES (?)",
        [(h,) for h in df["record_hash"]]
    )

    new_hashes = {
        row[0]
        for row in cursor.execute(f"""
            SELECT t.record_hash
            FROM {temp_table} t
            LEFT JOIN "{table_name}" m
                ON t.record_hash = m.record_hash
            WHERE m.record_hash IS NULL
        """)
    }

    cursor.execute(f"DROP TABLE {temp_table}")

    if not new_hashes:
        return 0

    new_df = df[df["record_hash"].isin(new_hashes)]

    # ------------------------------------------------------------
    # Insert new rows
    # ------------------------------------------------------------
    new_df = new_df.drop(columns=["exists"], errors="ignore")

    cols = list(new_df.columns)
    placeholders = ",".join("?" for _ in cols)
    col_str = ",".join(f'"{c}"' for c in cols)

    insert_sql = f"""
        INSERT INTO "{table_name}" ({col_str})
        VALUES ({placeholders})
    """
    rows = list(new_df.itertuples(index=False, name=None))

    try:
        conn.executemany(insert_sql, rows)
        conn.commit()

    except sql.IntegrityError as e:
        print("\n⚠️ Bulk insert failed, falling back to row-by-row")
        print("Reason:", e)

        inserted = 0
        for row in rows:
            try:
                conn.execute(insert_sql, row)
                inserted += 1
            except sql.IntegrityError as row_err:
                print("\n❌ Failed row:")
                print(dict(zip(new_df.columns, row)))
                print("➡️ Error:", row_err)


    # conn.executemany(
    #     insert_sql,
    #     list(new_df.itertuples(index=False, name=None))
    # )
    conn.commit()

    return len(new_df)

def vacuum_sqlite_db(db_path):
        """Vacuum (compact) the SQLite database to reduce file size and optimize performance."""
        if not os.path.exists(db_path):
            print(f"Database file '{db_path}' does not exist.")
            return
        conn = sql.connect(db_path)
        try:
            print(f"Vacuuming database: {db_path}")
            conn.execute("VACUUM")
            conn.commit()
            print("Vacuum completed successfully.")
        finally:
            conn.close()
            

def collect_mlst_data(input_dir):

    """Collect MLST data from TSV files in the input directory."""

    mlst_dir = os.path.join(input_dir, "mlst")
    if not os.path.exists(mlst_dir):
        print(f"Error: MLST directory '{mlst_dir}' does not exist.")
        return {}

    mlst_tsv_files = glob(os.path.join(mlst_dir, "*.tsv"))
    if not mlst_tsv_files:
        print(f"Error: No MLST TSV files found in '{mlst_dir}'.")
        return {}

    data = []
    for filename in mlst_tsv_files:
        sample_id = os.path.basename(filename).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)  # Use ID_MAPPER if available
        scheme, st, profile = None, None, None
        with open(filename, 'r') as f:
            mlst_line = f.readline()       
            parts = mlst_line.strip().split('\t')
            if len(parts) >= 3:
                scheme = parts[1] if parts[1] != '-' else None
                st = parts[2] if parts[2] != '-' else None
                # Get profile data (all fields after the ST)
                if len(parts) > 3:
                    profile = ','.join(parts[3:])
        data.append({
            'sample_id': sample_id,
            'scheme': scheme,
            'st': st,
            'profile': profile
        })

    df = pd.DataFrame(data)
    return df

def fetch_bracken_data(input_dir):
    """
    
    """
    
    bracken_dir = os.path.join(input_dir, "kraken2_bracken")
    if not os.path.exists(bracken_dir):
        print(f"Error: Bracken directory '{bracken_dir}' does not exist.")
        return {}
    bracken_tsv_files = glob(os.path.join(bracken_dir, "*.bracken.tsv"))
  
   
    bracken_data = []
    for filename in bracken_tsv_files:
        sample_id = os.path.basename(filename).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)  # Use ID_MAPPER if available
        df = pd.read_csv(filename, sep="\t", header=0)
       
        # Check if the dataframe is empty
        if df.empty:
            print(f"Warning: The file '{filename}' is empty. Skipping.")
            continue
        # Sort the dataframe by 'fraction_total_reads' in descending order and fetch the first row
        top_row = df.sort_values(by='fraction_total_reads', ascending=False).iloc[0]
        top_row_dict = top_row.to_dict()
        top_row_dict['sample_id'] = sample_id
        bracken_data.append(top_row_dict)
    

    df = pd.DataFrame(bracken_data)
    return df

def fetch_plasmidfinder_data(input_dir,conn):
    """
    def fetch_plasmidfinder_data(input_dir, conn):
        Fetches and processes PlasmidFinder data from the specified input directory and saves it to a database.
        This function searches for TSV files within the "plasmidfinder" subdirectory of the given input directory,
        processes the data, and inserts it into a database table named 'plasmidfinder'. If the table does not exist,
        it is created based on the structure of the first non-empty TSV file.
        Args:
            input_dir (str): The path to the input directory containing the "plasmidfinder" subdirectory.
            conn (sqlite3.Connection): A connection object to the SQLite database.
        Returns:
            dict: An empty dictionary is returned if the "plasmidfinder" directory does not exist.
                  Otherwise, the function does not return any value.
        Notes:
            - The function expects the TSV files to have specific column names, which are renamed to match
              the expected database schema.
            - The 'Query / Template length' column is split into two separate columns: 'qlen' and 'slen'.
            - Empty TSV files are skipped with a warning message.
            - If the "plasmidfinder" directory does not exist, an error message is printed, and the function exits.
    """

    plfinder_dir = os.path.join(input_dir, "plasmidfinder")
    if not os.path.exists(plfinder_dir):
        print(f"Error: The PlasmidFinder directory '{plfinder_dir}' does not exist.")
        return {}
    
    tsv_files = glob(os.path.join(plfinder_dir, "**", "*.tsv"), recursive=True)

    check_table = False
    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        # Rename columns to match the expected format
        df = df.rename(columns={
            'Database': 'database',
            'Plasmid': 'plasmid',
            'Identity': 'idty',
            'Query / Template length': 'qlen_slen',
            'Contig': 'contig',
            "Position in contig": "position",
            "Note": "note",
            "Accession number": "accession"
        })
        df['sample_id'] = sample_id
        # Split 'qlen_slen' into two separate columns and remove spaces
        df[['qlen', 'slen']] = df['qlen_slen'].str.replace(' ', '').str.split('/', expand=True)  # Remove unwanted columns
        df.drop(columns=['qlen_slen'], inplace=True)
        if not check_table:
            check_table = True
            create_table_if_not_exists(conn, 'plasmidfinder', df)  # Create table if it doesn't exist

        insert_data(conn, 'plasmidfinder', df)

def fetch_kleborate_data(input_dir):
    kleborate_dir = os.path.join(input_dir, "kleborate")
    if not os.path.exists(kleborate_dir):
        print(f"Error: The Kleborate directory '{kleborate_dir}' does not exist.")
        return pd.DataFrame()

    tsv_files = glob(os.path.join(kleborate_dir, "*.tsv"))
    if not tsv_files:
        print(f"Error: No TSV files found in '{kleborate_dir}'.")
        return pd.DataFrame()

    kleborate_data = []
    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        # Rename columns to match the expected format
        df = df.rename(columns={'strain': "sample_id"})
        df['sample_id'] = sample_id
        kleborate_data.append(df)

    if kleborate_data:
        return pd.concat(kleborate_data, ignore_index=True)
    else:
        return pd.DataFrame()

def fetch_vfdb_data(input_dir):
    vfdb_dir = os.path.join(input_dir, "diamond/VFDB")
    if not os.path.exists(vfdb_dir):
        print(f"Error: VFDB directory '{vfdb_dir}' does not exist.")
        return pd.DataFrame()
    
    tsv_files = glob(os.path.join(vfdb_dir, "*.tsv"))
    
    vfdb_data = []
    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        df["scov"] = 100. * round(df["length"] / df["slen"], 2)
        df["qcov"] = 100. * round(3 * df["length"] / df["qlen"], 2)
        df = df[(df["pident"] > 70) & (df["scov"] > 50) & (df["qcov"] > 50)]
        df["gene_name"] = df["stitle"].str.extract(r'\)\s*\((.*?)\)', expand=False)
        df["note"] = df["stitle"].str.replace(r'^.*?\)\s*\((.*?)\)\s*', '', regex=True)
        df["sample_id"] = sample_id

        vfdb_data.append(df)

    if vfdb_data:
        return pd.concat(vfdb_data, ignore_index=True)
    else:
        return pd.DataFrame()

def fetch_refbx_data(input_dir):
    refbx_dir = os.path.join(input_dir, "diamond")
    if not os.path.exists(refbx_dir):
        print(f"Error: diamond directory '{refbx_dir}' does not exist.")
        return pd.DataFrame()
    
    tsv_files = glob(os.path.join(refbx_dir, "**", "*__ref.diamond.tsv"), recursive=True)    
    refbx_data = []
    for tsv_file in tsv_files:
        
        sample_id,refacc,_,_ = os.path.basename(tsv_file).split('.')
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        refacc = refacc.split('__')[0].replace('v','.')

        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        df["scov"] = 100. * round(df["length"] / df["slen"], 2)
        df["qcov"] = 100. * round(3 * df["length"] / df["qlen"], 2)
        df = df[(df["pident"] > 60) ]       
        df["gene_name"] = df["stitle"].str.extract(r'\[gene=([^\]]+)\]', expand=False)
        df["locus_tag"] = df["stitle"].str.extract(r'\[locus_tag=([^\]]+)\]', expand=False)
        df["product"] = df["stitle"].str.extract(r'\[protein=([^\]]+)\]', expand=False)
        df["sample_id"] = sample_id
        df["ref_acc"] = refacc
        df.drop(columns=['stitle'], inplace=True, errors='ignore')
        refbx_data.append(df)

    if refbx_data:
        return pd.concat(refbx_data, ignore_index=True)
    else:
        return pd.DataFrame()

def fetch_rgi_data(input_dir):
    rgi_dir = os.path.join(input_dir, "rgi")
    if not os.path.exists(rgi_dir):
        print(f"Error: The RGI directory '{rgi_dir}' does not exist.")
        return pd.DataFrame

    tsv_files = glob(os.path.join(rgi_dir, "**", "*.rgi.txt"), recursive=True)
    if not tsv_files:
        print(f"Error: No RGI TSV files found in '{rgi_dir}'.")
        return pd.DataFrame()

    rgi_data = []

    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        
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

        df["sample_id"] = sample_id
        rgi_data.append(df)
    
    if rgi_data:
        return pd.concat(rgi_data, ignore_index=True)
    else:
        return pd.DataFrame()

def fetch_amrfinder_data(input_dir):

    amrfinder_dir = os.path.join(input_dir, "amrfinder")
    if not os.path.exists(amrfinder_dir):
        print(f"Error: The AMRFinder directory '{amrfinder_dir}' does not exist.")
        return pd.DataFrame()
    tsv_files = []
    tsv_files = glob(os.path.join(amrfinder_dir, "*.amrfinder.tsv"))
    tsv_files += glob(os.path.join(amrfinder_dir, "*.mutations.tsv"))

    if not tsv_files:
        print(f"Error: No AMRFinder TSV files found in '{amrfinder_dir}'.")
        return pd.DataFrame()

    amrfinder_data = []
    for tsv_file in tsv_files:
        sample_id = os.path.basename(tsv_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)

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
            "Protein identifier", "HMM id", "HMM description"
        ]
        df = df.drop(columns=columns_to_drop, errors='ignore')

        # Rename columns for consistency
        column_renames = {
            "Contig id": "contig_id",
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
            "Accession of closest sequence": "accession_of_closest_sequence",
            "Name of closest sequence": "name_of_closest_sequence",
        }
        df = df.rename(columns=column_renames)
        # Add sample_id column
        df["sample_id"] = sample_id
        amrfinder_data.append(df)

    return pd.concat(amrfinder_data, ignore_index=True).drop_duplicates() if amrfinder_data else pd.DataFrame()

def fetch_snippy_data(input_dir):

    snippy_dir = os.path.join(input_dir, "snippy")
    if not os.path.exists(snippy_dir):
        print(f"Error: The snippy directory '{snippy_dir}' does not exist.")
        return pd.DataFrame()
    
    csv_files = glob(os.path.join(snippy_dir, "**", "*.csv"), recursive=True)
    if not csv_files:
        print(f"Error: No snippy CSV files found in '{snippy_dir}'.")
        return pd.DataFrame()
    
    snippy_data = []
    for csv_file in csv_files:
        sample_id = os.path.basename(csv_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)

        try:
            df = pd.read_csv(csv_file, header=0)
        except Exception as e:
            print(f"Error reading file '{csv_file}': {e}")
            continue

        if df.empty:
            print(f"Warning: The file '{csv_file}' is empty. Skipping.")
            continue

        if "EFFECT" in df.columns:
            df = df[~df["EFFECT"].str.contains("synonymous_variant", na=False)]

        df["sample_id"] = sample_id
        snippy_data.append(df)

    return pd.concat(snippy_data, ignore_index=True) if snippy_data else pd.DataFrame()

def fetch_snippy_genome(input_dir,max_workers=1):
    snippy_dir = os.path.join(input_dir, "snippy_contig")
    if not os.path.exists(snippy_dir):
        print(f"Error: The snippy directory '{snippy_dir}' does not exist.")
        return pd.DataFrame()
    
    vcf_files = glob(os.path.join(snippy_dir, "**", "*.snpEff.vcf"), recursive=True)
    if not vcf_files:
        print(f"Error: No snippy.snpEff vcf files found in '{snippy_dir}'.")
        return pd.DataFrame()
    
    df_lst = process_vcf_files_parallel(vcf_files,max_workers=max_workers)
    if df_lst:
        df = pd.concat(df_lst, ignore_index=True)
        if "ANN_Annotation" in df.columns:
            df = df[~df['ANN_Annotation'].str.contains("synonymous_variant", na=False)]
        if "ANN_Gene_ID" in df.columns:
            df = df[df['ANN_Gene_ID'].astype(str).str.len() > 1]
        
        if "sample_id" in df.columns:
            df["sample_id"] = df["sample_id"].apply(lambda x: ID_MAPPER.get(x, x))

        return df
    else:
        return pd.DataFrame()

def fetch_quast_data(input_dir):

    quast_dir = os.path.join(input_dir, "quast")
    if not os.path.exists(quast_dir):
        print(f"Error: The QUAST directory '{quast_dir}' does not exist.")
        return pd.DataFrame()
    
    tsv_files = glob(os.path.join(quast_dir, "**", "transposed_report.tsv"), recursive=True)
    if not tsv_files:
        print(f"Error: No QUAST TSV files found in '{quast_dir}'.")
        return pd.DataFrame()
    quast_data = []
    for tsv_file in tsv_files:
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue

        if "Assembly" in df.columns:
            sample_id = str(df["Assembly"].iloc[0]).split('.')[0]
        else:
            sample_id = os.path.basename(os.path.dirname(tsv_file))

        sample_id = ID_MAPPER.get(sample_id, sample_id)

        # Rename columns for consistency
        column_renames = {
            "Assembly": "assembly",
            "# contigs (>= 0 bp)": "contigs_0bp",
            "# contigs (>= 1000 bp)": "contigs_1kbp",
            "# contigs (>= 5000 bp)": "contigs_5kbp",
            "# contigs (>= 10000 bp)": "contigs_10kbp",
            "# contigs (>= 25000 bp)": "contigs_25kbp",
            "# contigs (>= 50000 bp)": "contigs_50kbp",
            "Total length (>= 0 bp)": "total_length_0bp",
            "Total length (>= 1000 bp)": "total_length_1kbp",
            "Total length (>= 5000 bp)": "total_length_5kbp",
            "Total length (>= 10000 bp)": "total_length_10kbp",
            "Total length (>= 25000 bp)": "total_length_25kbp",
            "Total length (>= 50000 bp)": "total_length_50kbp",
            "# contigs": "contigs",
            "Largest contig": "largest_contig",
            "Total length": "total_length",
            "Reference length": "reference_length",
            "GC (%)": "gc_percent",
            "Reference GC (%)": "reference_gc_percent",
            "N50": "n50",
            "NG50": "ng50",
            "N90": "n90",
            "NG90": "ng90",
            "auN": "aun",
            "auNG": "aung",
            "L50": "l50",
            "LG50": "lg50",
            "L90": "l90",
            "LG90": "lg90",
            "# misassemblies": "misassemblies",
            "# misassembled contigs": "misassembled_contigs",
            "Misassembled contigs length": "misassembled_contigs_length",
            "# local misassemblies": "local_misassemblies",
            "# scaffold gap ext. mis.": "scaffold_gap_ext_mis",
            "# scaffold gap loc. mis.": "scaffold_gap_loc_mis",
            "# unaligned mis. contigs": "unaligned_mis_contigs",
            "# unaligned contigs": "unaligned_contigs",
            "Unaligned length": "unaligned_length",
            "Genome fraction (%)": "genome_fraction_percent",
            "Duplication ratio": "duplication_ratio",
            "# N's per 100 kbp": "ns_per_100kbp",
            "# mismatches per 100 kbp": "mismatches_per_100kbp",
            "# indels per 100 kbp": "indels_per_100kbp",
            "Largest alignment": "largest_alignment",
            "Total aligned length": "total_aligned_length",
            "NA50": "na50",
            "NGA50": "nga50",
            "NA90": "na90",
            "NGA90": "nga90",
            "auNA": "auna",
            "auNGA": "aunga",
            "LA50": "la50",
            "LGA50": "lga50",
            "LA90": "la90",
            "LGA90": "lga90"
        }
        df.rename(columns=column_renames, inplace=True)
        df["sample_id"] = sample_id
        
        quast_data.append(df)

    return pd.concat(quast_data, ignore_index=True) if quast_data else pd.DataFrame()

def fetch_orfs_seq(input_dir,output_dir, db_name):
    prokka_dir = os.path.join(input_dir, "prokka")
    if not os.path.exists(prokka_dir):
        print(f"Error: The Prokka directory '{prokka_dir}' does not exist.")
        return pd.DataFrame()

    # Connect to SQLite database (or create it)
    name,ext = db_name.split('.')
    db_path = os.path.join(output_dir, f"{name}_sequence.{ext}")
    conn = sql.connect(db_path)
    ffn_files = glob(os.path.join(prokka_dir, "**", "*.ffn"), recursive=True)
    check_db = True
    for ffn in ffn_files:
        sample_id = os.path.basename(ffn).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        print(f'\t">>Fetching {sample_id} Open Reading Frames.')
        orf_data = {}
        with open(ffn) as hdl:
            for rec in SeqIO.parse(hdl,'fasta'):
                orf_data[rec.id] = {"sample_id":sample_id,"locus_tag":rec.id, 
                                    'dscr': str(rec.description), 'nt_seq':str(rec.seq), 'aa_seq':None}
        
        faa = os.path.join(os.path.dirname(ffn),f'{sample_id}.faa')
        with open(faa) as hdl:
            for rec in SeqIO.parse(hdl,'fasta'):
                if rec.id not in orf_data:
                    orf_data[rec.id] = {"sample_id":sample_id, "locus_tag":rec.id, 
                                    'dscr': str(rec.description), 'nt_seq':None, 'aa_seq':None}
                orf_data[rec.id]['aa_seq'] = str(rec.seq)

        orf_df = pd.DataFrame(orf_data.values())
        if check_db:
            check_db = False
            create_table_if_not_exists(conn, 'orfs', orf_df)

        insert_data(conn, 'orfs', orf_df)
       
    
    conn.close()
    return db_path

def fetch_prokka_data(input_dir):

    prokka_dir = os.path.join(input_dir, "prokka")
    if not os.path.exists(prokka_dir):
        print(f"Error: The Prokka directory '{prokka_dir}' does not exist.")
        return pd.DataFrame()

    tsv_files = glob(os.path.join(prokka_dir, "**", "*.tsv"), recursive=True)
    if not tsv_files:
        print(f"Error: No Prokka TSV files found in '{prokka_dir}'.")
        return pd.DataFrame()

    prokka_data = []
    for tsv_file in tsv_files:
        print(f"\t>>{tsv_file}")
        sample_id = os.path.basename(tsv_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        df = pd.read_csv(tsv_file, sep="\t", header=0)
        if df.empty:
            print(f"Warning: The file '{tsv_file}' is empty. Skipping.")
            continue
        df["sample_id"] = sample_id
        prokka_data.append(df)
    
    return pd.concat(prokka_data, ignore_index=True) if prokka_data else pd.DataFrame() 

def fetch_gene_diff(input_dir):
    genediff_dir = os.path.join(input_dir, "genediff")
    if not os.path.exists(genediff_dir):
        print(f"Error: The GeneDiff directory '{genediff_dir}' does not exist.")
        return pd.DataFrame()

    json_files = glob(os.path.join(genediff_dir, "*.json"))
    if not json_files:
        print(f"Error: No GeneDiff JSON files found in '{genediff_dir}'.")
        return pd.DataFrame()

    genediff_data = []
    for json_file in json_files:
        sample_id = os.path.basename(json_file).split('.')[0]
        sample_id = ID_MAPPER.get(sample_id, sample_id)
        with open(json_file, 'r') as f:
            gene_list =  json.load(f)
            for gene in gene_list:
                
                info = gene.get('gene_info',None)    
                qry_id = gene.get('query_sequence_id',None)
                
                if info is None or qry_id is None:
                    print(f"Warning: Missing gene_info or query_sequence_id in {json_file}. Skipping this gene.")
                    continue

                effect_types = None
                genediff_data.append({
                    "sample_id": sample_id,
                    "qry_id": qry_id,
                    "gene_name": info.get("gene"),
                    "locus_tag": info.get("locus_tag"),
                    "product": info.get("protein"),
                    "accession": info.get("protein_id"),
                    "identity": gene.get("blast_identity"),
                    "scov": gene.get("blast_scoverage"),
                    "qcov": gene.get("blast_qcoverage"),
                    "num_mutations": gene.get("num_mutations"),
                    "effects": effect_types,
                    "mutations": None,
                    "ref_seq": None,
                    "qry_seq": None,
                })

                # Add mutations if available
                protein_analysis = gene.get("protein_analysis", None)
                if protein_analysis:
                    protein_changes = protein_analysis.get("protein_changes", [])
                    genediff_data[-1]["mutations"] = '; '.join(protein_changes)
                    effect_types = protein_analysis.get("effect_types", None)
                    genediff_data[-1]["effects"] = '; '.join(effect_types)

                    # Only add ref/qry_seq if effect is not exclusively synonymous or missense
                    effect_types = set(effect_types) if effect_types else set()
                    effect_types.difference_update({"synonymous", "missense"})
                    if len(effect_types) > 0:
                        genediff_data[-1]["ref_seq"] = protein_analysis.get("ref_protein")
                        genediff_data[-1]["qry_seq"] = protein_analysis.get("query_protein")
    
    return pd.DataFrame(genediff_data) if genediff_data else pd.DataFrame()

def save_to_db(input_dir, output_dir, db_name, add_seq = False,cpus=1,skip=[]):
    """Main function to save data to SQLite database."""
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Connect to SQLite database (or create it)
    db_path = os.path.join(output_dir, db_name)
    conn = sql.connect(db_path)

    print("\n=== Processing MLST Data ===")
    mlst_df = collect_mlst_data(input_dir)
    if not mlst_df.empty:
        create_table_if_not_exists(conn, 'mlst', mlst_df)
        inserted_mlst = insert_data(conn, 'mlst', mlst_df)
        print(f"MLST data processed: {len(mlst_df)} records found, {inserted_mlst} new records inserted.")
    else:
        print("No MLST data found to process.")

    print("\n=== Processing Bracken Data ===")
    bracken_df = fetch_bracken_data(input_dir)
    if not bracken_df.empty:
        create_table_if_not_exists(conn, 'bracken', bracken_df)
        inserted_bracken = insert_data(conn, 'bracken', bracken_df)
        print(f"Bracken data processed: {len(bracken_df)} records found, {inserted_bracken} new records inserted.")
    else:
        print("No Bracken data found to process.")

    print("\n=== Processing PlasmidFinder Data ===")
    fetch_plasmidfinder_data(input_dir,conn)
    print("PlasmidFinder data processed.")

    print("\n=== Processing Kleborate Data ===")
    kleborate_df = fetch_kleborate_data(input_dir)
    create_table_if_not_exists(conn, 'kleborate', kleborate_df)
    insert_data(conn, 'kleborate', kleborate_df)
    print(f"Kleborate data processed: {len(kleborate_df)} records found.")

    print("\n=== Processing VFDB Data ===")
    vfdb_df = fetch_vfdb_data(input_dir)
    create_table_if_not_exists(conn, 'vfdb', vfdb_df)
    insert_data(conn, 'vfdb', vfdb_df)
    print(f"VFDB data processed: {len(vfdb_df)} records found.")
    
    print("\n=== Processing Ref Blastx Data ===")
    refbx_df = fetch_refbx_data(input_dir)
    create_table_if_not_exists(conn, 'refbx', refbx_df)
    insert_data(conn, 'refbx', refbx_df)
    print(f"RefBx data processed: {len(refbx_df)} records found.")
    
    print("\n=== Processing RGI Data ===")
    rgi_df = fetch_rgi_data(input_dir)
    create_table_if_not_exists(conn, 'rgi', rgi_df)
    insert_data(conn, 'rgi', rgi_df)
    print(f"RGI data processed: {len(rgi_df)} records found.")

    print("\n=== Processing AMRFinder Data ===")
    amrfinder_df = fetch_amrfinder_data(input_dir)
    create_table_if_not_exists(conn, 'amrfinder', amrfinder_df)
    insert_data(conn, 'amrfinder', amrfinder_df)
    print(f"AMRFinder data processed: {len(amrfinder_df)} records found.")

    print("\n=== Processing QUAST Data ===")
    quast_df = fetch_quast_data(input_dir)
    create_table_if_not_exists(conn, 'quast', quast_df)
    insert_data(conn, 'quast', quast_df)
    print(f"QUAST data processed: {len(quast_df)} records found.")

    print("\n=== Processing Prokka Data ===")
    prokka_df = fetch_prokka_data(input_dir)
    create_table_if_not_exists(conn, 'prokka', prokka_df)
    insert_data(conn, 'prokka', prokka_df)
    print(f"Prokka data processed: {len(prokka_df)} records found.")
    
    if('snippy' not in skip):
        print("\n=== Processing Snippy Data ===")
        snippy_df = fetch_snippy_data(input_dir)
        create_table_if_not_exists(conn, 'snippy', snippy_df)
        insert_data(conn, 'snippy', snippy_df,add_timestamp=False)
        print(f"Snippy data processed: {len(snippy_df)} records found.")
    else:
        print("\n>>> Skipping Snippy Data Processing <<<")
            
    if('snippy_genome' not in skip):
        print("\n=== Processing Snippy Genome Data ===")
        snippy_genome_df = fetch_snippy_genome(input_dir,cpus)
        create_table_if_not_exists(conn, 'snippy_genome', snippy_genome_df)
        insert_data(conn, 'snippy_genome', snippy_genome_df, add_timestamp=False)
        print(f"Snippy genome data processed: {len(snippy_genome_df)} records found.")
    else:
        print("\n>>> Skipping Snippy Genome Data Processing <<<")

    if('genediff' not in skip):
        print("\n=== Processing GeneDiff Data ===")
        genediff_df = fetch_gene_diff(input_dir)
        create_table_if_not_exists(conn, 'genediff', genediff_df)
        insert_data(conn, 'genediff', genediff_df, add_timestamp=False)
        print(f"GeneDiff data processed: {len(genediff_df)} records found.")
    else:
        print("\n>>> Skipping GeneDiff Data Processing <<<")

    if add_seq:
        print("\n=== Processing ORFs Sequences ===")
        seqdb = fetch_orfs_seq(input_dir,output_dir,db_name)
        print(f"ORFs sequences were saved in: {seqdb}")

    # Close the database connection
    conn.close()
    print("Data saved to database successfully.")

    print("\n*** vacuume database ***")
    vacuum_sqlite_db(db_path)
 
def load_meta_data(meta_name, meta_data,output_dir, db_name):
    """
    Loads metadata from a CSV file into a specified table in the SQLite database.
    The CSV file should have a header row with column names matching the table schema.
    """
    db_path = os.path.join(output_dir, db_name)
    conn = sql.connect(db_path)
    
    df = pd.read_csv(meta_data,header=0, dtype=str)
    df = df.fillna('')  # Fill NaN values with empty strings for consistency

    if df.empty:
        print(f"Warning: The metadata file '{meta_data}' is empty. No data loaded.")
        return
    df = df.drop_duplicates()
    # Create table if it doesn't exist
    create_table_if_not_exists(conn, meta_name, df)

    # Insert data into the table
    inserted_rows = insert_data(conn, meta_name, df)
    print(f"Metadata loaded into '{meta_name}' table: {inserted_rows} rows inserted.")

    conn.close()

def read_csv_to_dict(csv_file):
    """
    Reads a CSV file with two columns and returns a dictionary.
    The first column is used as the key, the second as the value.
    """
    result = {}
    with open(csv_file, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) >= 2:
                key, value = parts[0], parts[1]
                result[key] = value
    return result

def main():
    parser = argparse.ArgumentParser(description="Save data to database.")
    parser.add_argument('--input_dir', required=True, help='Path to the input directory')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory')
    parser.add_argument('--db_name', required=False, default='trackpath_results.db', help='Name of the SQLite database file')
    parser.add_argument('--add_seq', action='store_true', help='If provided, perform additional sequence-related processing')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use for parallel processing')
    parser.add_argument('--id_mapper', type=str, default=None, help='CSV file mapping old_id to new_id (columns: old_id,new_id) without header')
    parser.add_argument('--meta_name', type=str, default=None, help='Name of the metadata table to create')
    parser.add_argument('--meta_data', type=str, default=None, help='CSV file containing metadata for the table')
    parser.add_argument('--skip', type=str, default='snippy_genome',
                        help="Comma-separated list of results to skip saving (snippy, snippy_genome, genediff). Example: snippy,snippy_genome")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir
    db_name = args.db_name
   
    # Check that if meta_name or meta_data is used, both should be provided
    
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Database name: {db_name}")
    print(f"The database contains ORF sequences? {['No','Yes'][args.add_seq]}" )
    print(f"skipping results: {args.skip}")

    if args.id_mapper:
        global ID_MAPPER
        
        ID_MAPPER = read_csv_to_dict(args.id_mapper)
        print(f"ID mapper loaded with {len(ID_MAPPER)} mappings.")

    if (args.meta_name is not None) ^ (args.meta_data is not None):
        parser.error("Both --meta_name and --meta_data must be provided together.")
    
    if args.meta_name and args.meta_data:
        print(f"Loading metadata into table '{args.meta_name}' from file '{args.meta_data}'")
        load_meta_data(args.meta_name, args.meta_data, output_dir, db_name)
        

    save_to_db(input_dir, output_dir, db_name,args.add_seq,args.cpus,skip=args.skip.split(','))
    print("Data processing complete.")

if __name__ == "__main__":
    main()