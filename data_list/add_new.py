# --*-- conding:utf-8 --*--
# @Time : 3/26/25 8:48 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : add_new.py

# Mapping dictionary from three-letter residue codes to one-letter codes.
THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def convert_three_to_one(seq):
    """
    Convert a residue sequence from hyphen-separated three-letter codes to a single-letter code string.

    Parameters:
        seq (str): e.g., "ILE-LEU-MET-GLU-LEU-MET-ALA-GLY-GLY-ASP-LEU-LYS-SER-PHE"

    Returns:
        str: e.g., "ILMELMAGGDLKSF"
    """
    parts = seq.strip().split('-')
    one_letter = ''
    for p in parts:
        code = p.strip().upper()
        one_letter += THREE_TO_ONE.get(code, 'X')  # Use 'X' for unknown residues
    return one_letter


def merge_and_convert(original_path, new_path, output_path=None):
    """
    Merge the contents of the original file and new file.
    For each line in the new file, convert the residue sequence (last column)
    from three-letter codes to a one-letter code string.
    The merged content is written to output_path. If output_path is None,
    the original file is overwritten.

    Parameters:
        original_path (str): Path to the original file.
        new_path (str): Path to the new file.
        output_path (str, optional): Path for the merged output file.
    """
    # Read the original file
    with open(original_path, 'r', encoding='utf-8') as f:
        original_lines = f.readlines()

    # Process the new file: convert the last column (residue sequence)
    new_processed_lines = []
    with open(new_path, 'r', encoding='utf-8') as f:
        for line in f:
            # Remove trailing newline
            line = line.rstrip("\n")
            # Assume columns are separated by whitespace or tabs.
            # Here we split by tab first. Adjust if needed.
            columns = line.split("\t")
            if len(columns) >= 5:
                # The last column is the residue sequence in three-letter codes.
                original_seq = columns[-1]
                converted_seq = convert_three_to_one(original_seq)
                # Replace the last column with the converted sequence.
                columns[-1] = converted_seq
                # Rejoin with tab (you can also use a specific delimiter if needed)
                new_line = "\t".join(columns)
            else:
                new_line = line  # Leave the line unchanged if it doesn't have enough columns
            new_processed_lines.append(new_line + "\n")

    # Determine output file path: if not provided, overwrite the original file.
    if output_path is None:
        output_path = original_path

    # Write the merged content: first original file's lines, then the processed new lines.
    with open(output_path, 'w', encoding='utf-8') as f:
        for line in original_lines:
            f.write(line)
        # Optionally, insert a blank line between sections:
        # f.write("\n")
        for line in new_processed_lines:
            f.write(line)

    print(f"Merged file written to: {output_path}")


if __name__ == "__main__":
    # Set the paths to your original and new files
    original_file = "integrated_data.txt"  # Path to the original text file
    new_file = "14_14.txt"  # Path to the new text file to be appended
    output_file = "all_data_list.txt"  # Output file path (set to original_file to overwrite)

    merge_and_convert(original_file, new_file, output_file)
