# --*-- conding:utf-8 --*--
# @Time : 3/23/25 11:15 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : sum_list.py

if __name__ == '__main__':

    amino_acid_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    # List of input files based on chain length categories
    input_files = ['5_7.txt', '8_9.txt', '10_11.txt', '12_13.txt']

    # Name of the output file
    output_file = 'integrated_data.txt'

    # Open the output file in write mode with UTF-8 encoding
    with open(output_file, 'w', encoding='utf-8') as out_f:
        # Process each input file
        for file_name in input_files:
            try:
                # Open the input file in read mode with UTF-8 encoding
                with open(file_name, 'r', encoding='utf-8') as in_f:
                    for line in in_f:
                        # Remove leading and trailing whitespace
                        line = line.strip()
                        # Split the line into fields using tab as separator
                        fields = line.split('\t')
                        # Check if the line has exactly 5 fields
                        if len(fields) != 5:
                            print(f"Warning: Line in {file_name} does not have 5 fields: {line}")
                            continue
                        # Extract the amino acid sequence (last field)
                        sequence = fields[4]
                        # Split the sequence by hyphens and clean each code
                        three_letter_codes = [code.strip().upper() for code in sequence.split('-')]
                        # Convert to single-letter codes
                        single_letter_sequence = ''
                        for code in three_letter_codes:
                            if code in amino_acid_map:
                                single_letter_sequence += amino_acid_map[code]
                            else:
                                single_letter_sequence += '?'
                                print(f"Warning: Unknown amino acid code {code} in {file_name}")
                        # Replace the original sequence with the single-letter sequence
                        fields[4] = single_letter_sequence
                        # Write the modified line to the output file
                        out_f.write('\t'.join(fields) + '\n')
            except FileNotFoundError:
                print(f"Warning: File {file_name} not found.")