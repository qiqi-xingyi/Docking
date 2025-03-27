# --*-- conding:utf-8 --*--
# @Time : 3/26/25 8:34 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : chain_14_pre.py

import os
from docking import Fileprepare

grouped_result_dir = "grouped_result/chain_14"
selected_dir = "selected"
output_protein_dir = "14_protein_pdbqt"
output_ligand_dir = "14_ligand_pdbqt"

os.makedirs(output_protein_dir, exist_ok=True)
os.makedirs(output_ligand_dir, exist_ok=True)

def extract_protein_id(xyz_filename):
    base = os.path.splitext(xyz_filename)[0]
    parts = base.split('_')
    return parts[0]

if __name__ == '__main__':

    for protein_dir in os.listdir(grouped_result_dir):
        protein_path = os.path.join(grouped_result_dir, protein_dir)
        if not os.path.isdir(protein_path):
            continue

        for file_name in os.listdir(protein_path):
            if not file_name.endswith(".xyz"):
                continue

            if "_top_" not in file_name:
                protein_id_temp = extract_protein_id(file_name)
                top_1_filename = f"{protein_id_temp}_top_1.xyz"
                top_1_path = os.path.join(protein_path, top_1_filename)
                if os.path.exists(top_1_path):
                    print(f"Skipping {file_name} because {top_1_filename} exists.")
                    continue

            xyz_path = os.path.join(protein_path, file_name)
            protein_id = extract_protein_id(file_name)


            protein_output_dir = os.path.join(output_protein_dir, protein_id)
            os.makedirs(protein_output_dir, exist_ok=True)

            print(f"\n=== Processing protein '{protein_id}' from '{file_name}' ===")
            docking_obj = Fileprepare(
                xyz_file=xyz_path,
                docking_folder=protein_output_dir
            )
            docking_obj.run_pipeline()


            ligand_dir = os.path.join(selected_dir, protein_id)
            ligand_mol2 = os.path.join(ligand_dir, f"{protein_id}_ligand.mol2")

            if os.path.exists(ligand_mol2):
                ligand_output_dir = os.path.join(output_ligand_dir, protein_id)
                os.makedirs(ligand_output_dir, exist_ok=True)

                print(f"=== Processing ligand for '{protein_id}' from '{ligand_mol2}' ===")
                docking_obj.prepare_ligand_pdbqt(
                    ligand_mol2_path=ligand_mol2,
                    output_folder=ligand_output_dir
                )
            else:
                print(f"Warning: ligand file not found for {protein_id}, expected at {ligand_mol2}")

    print("\nAll processing done!")
