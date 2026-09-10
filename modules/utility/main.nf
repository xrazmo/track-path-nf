

process EXTRACT_PROTEINS {
    tag "${genbank_file}"
    // Process to extract proteins from GenBank files
    input:
    tuple val(meta) , path(genbank_file)

    // Output definition
    output:
    path "*.faa", emit: proteins
    path "*.log", emit: log

    // Shell script for protein extraction
    script:
    // Generate output filename based on input
    def prefix = genbank_file.baseName
    """
    # Protein extraction script
    awk '
    # State tracking variables
    /^     CDS / {in_cds = 1; protein_line = ""; header = ""; locus_tag = ""; gene_name = ""; product = ""}
    /^     gene / {in_cds = 0}
    /^     misc_feature / {in_cds = 0}

    # Capture the gene name
    in_cds && /^                     \\/gene=/ {
        # Remove quotes and clean up
        gsub(/^[ \\t]*\\/gene="|"$/, "", \$0)
        gsub(/"/, "", \$0)  # Remove any remaining quotes
        gene_name = \$0
    }

    # Capture the locus_tag
    in_cds && /^                     \\/locus_tag=/ {
        # Remove quotes and clean up
        gsub(/^[ \\t]*\\/locus_tag="|"$/, "", \$0)
        gsub(/"/, "", \$0)  # Remove any remaining quotes
        locus_tag = \$0
    }

    # Capture the protein name (optional)
    in_cds && /^                     \\/product=/ {
        # Remove quotes and clean up
        gsub(/^[ \\t]*\\/product="|"$/, "", \$0)
        gsub(/"/, "", \$0)  # Remove any remaining quotes
        product = \$0
    }

    # Capture multi-line translation
    in_cds && /^                     \\/translation=/ {
        # Remove quotes and leading/trailing whitespace
        gsub(/^[ \\t]*\\/translation="|"$/, "", \$0)
        gsub(/"/, "", \$0)  # Remove any remaining quotes
        protein_line = \$0
        in_translation = 1
        next
    }

    # Continue capturing translation across multiple lines
    in_translation && /^                     [a-zA-Z]+/ {
        # Remove leading whitespace and continue building protein sequence
        gsub(/^[ \\t]*/, "", \$0)
        gsub(/"/, "", \$0)  # Remove any quotes in translation
        protein_line = protein_line \$0
        next
    }

    # End of translation block
    in_translation && !/^                     [a-zA-Z]+/ {
        in_translation = 0
        # Only output if we have a protein sequence
        if (protein_line != "") {
            # Construct header with multiple components
            header_parts[0] = locus_tag
            header_parts_count = 1
            
            # Add gene name if available
            if (gene_name != "") {
                header_parts[header_parts_count++] = gene_name
            }
            
            # Add product name if available
            if (product != "") {
                header_parts[header_parts_count++] = product
            }
            
            # Join header parts
            header = header_parts[0]
            for (i = 1; i < header_parts_count; i++) {
                header = header " | " header_parts[i]
            }
            
            # Default headers if no identifier found
            if (header == "") header = "Unnamed_Protein"
            
            print ">" header
            print protein_line
            
            # Reset variables
            protein_line = ""
            header = ""
            locus_tag = ""
            gene_name = ""
            product = ""
            delete header_parts
        }
    }

    END {
        # Capture the last protein if exists
        if (protein_line != "") {
            # Construct header with multiple components
            header_parts[0] = locus_tag
            header_parts_count = 1
            
            # Add gene name if available
            if (gene_name != "") {
                header_parts[header_parts_count++] = gene_name
            }
            
            # Add product name if available
            if (product != "") {
                header_parts[header_parts_count++] = product
            }
            
            # Join header parts
            header = header_parts[0]
            for (i = 1; i < header_parts_count; i++) {
                header = header " | " header_parts[i]
            }
            
            # Default headers if no identifier found
            if (header == "") header = "Unnamed_Protein"
            
            print ">" header
            print protein_line
        }
    }' ${genbank_file} | \
    awk '
    BEGIN { 
        line_length = 60  # Break protein sequences into 60-character lines
    }
    /^>/ { 
        if (NR > 1) print ""  # Add blank line between entries
        print \$0  # Print header
        next 
    }
    {
        # Break long sequences into specified line length
        while (length(\$0) > line_length) {
            print substr(\$0, 1, line_length)
            \$0 = substr(\$0, line_length + 1)
        }
        print \$0
    }
    END { print "" }  # Ensure file ends with newline
    ' > ${prefix}.faa 2> ${prefix}_protein_extraction.log


    # Check if any proteins were extracted
    if [ ! -s ${prefix}.faa ]; then
        echo "WARNING: No proteins extracted from ${genbank_file}" >> ${prefix}_protein_extraction.log
    fi
    """
}
