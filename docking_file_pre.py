# --*-- conding:utf-8 --*--
# @Time : 3/20/25 8:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking_file_pre.py

import os
from docking import Fileprepare

grouped_result_dir = "grouped_result"
selected_dir = "selected"
output_protein_dir = "protein_pdbqt"
output_ligand_dir = "ligand_pdbqt"

os.makedirs(output_protein_dir, exist_ok=True)
os.makedirs(output_ligand_dir, exist_ok=True)


def extract_protein_id(xyz_filename):
    base = os.path.splitext(xyz_filename)[0]
    parts = base.split('_')
    return parts[0]


if __name__ == '__main__':

    for chain_dir in os.listdir(grouped_result_dir):
        chain_path = os.path.join(grouped_result_dir, chain_dir)
        if not os.path.isdir(chain_path):
            continue

        for protein_dir in os.listdir(chain_path):
            protein_path = os.path.join(chain_path, protein_dir)
            if not os.path.isdir(protein_path):
                continue

            # e.g. grouped_result/chain_5/1fkf/*.xyz
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











