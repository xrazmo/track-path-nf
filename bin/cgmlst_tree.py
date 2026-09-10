import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
import argparse

# Parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate cgMLST tree from allelic distance matrix.")
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to distances.tsv file containing allelic distances."
    )
    parser.add_argument(
        "-o", "--output",
        default="cgmlst_tree.nwk",
        help="Path to output Newick file (default: cgmlst_tree.nwk)."
    )
    return parser.parse_args()

# Read the distance matrix from distances.tsv
def read_distance_matrix(file_path):
    # Read the tab-separated distance matrix
    df = pd.read_csv(file_path, sep='\t', index_col=0)
    # Ensure the matrix is symmetric and square
    if not df.shape[0] == df.shape[1]:
        raise ValueError("Distance matrix must be square")
    # Convert to numpy array for clustering
    distance_matrix = df.values
    # Ensure the matrix is symmetric
    if not np.allclose(distance_matrix, distance_matrix.T):
        raise ValueError("Distance matrix must be symmetric")
    return df, distance_matrix

# Perform hierarchical clustering (UPGMA)
def compute_linkage(distance_matrix):
    # Convert distance matrix to condensed form (upper triangle)
    condensed_dist = squareform(distance_matrix, checks=False)
    # Perform UPGMA clustering
    Z = linkage(condensed_dist, method='average')  # 'average' corresponds to UPGMA
    return Z

# Convert linkage to a Newick tree
def linkage_to_newick(Z, labels):
    # Convert linkage matrix to a tree
    tree = to_tree(Z, rd=False)
    
    # Recursive function to build Newick string
    def to_newick(node, labels):
        if node.is_leaf():
            return labels[node.id]
        left = to_newick(node.left, labels)
        right = to_newick(node.right, labels)
        d = round(node.dist/2,4)
        return f"({left}:{d},{right}:{d})"
    
    newick = to_newick(tree, labels) + ";"
    return newick

# Main function
def main():
    # Parse arguments
    args = parse_arguments()
    distances_file = args.input
    output_file = args.output
    
    # Read distance matrix
    df, distance_matrix = read_distance_matrix(distances_file)
    isolate_names = df.index.tolist()
    
    # Compute hierarchical clustering
    Z = compute_linkage(distance_matrix)
    
    # Convert to Newick format
    newick = linkage_to_newick(Z, isolate_names)
    
    # Save the Newick string to a file
    with open(output_file, "w") as f:
        f.write(newick)
    print(f"Newick tree saved as {output_file}")

if __name__ == "__main__":
    main()