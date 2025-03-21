# --*-- conding:utf-8 --*--
# @Time : 3/20/25 8:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking.py

import os
import glob

from docking import Docking

def batch_docking_for_xyz(chain_dir, ligand_mol2, output_base_dir):
    """
    Loop over all XYZ files in 'chain_dir', generate a full-atom PDB for each,
    and run AutoDock Vina docking with the given MOL2 ligand.

    Parameters:
    - chain_dir (str): Path to the directory containing .xyz files (e.g., 'grouped_result/chain_5')
    - ligand_mol2 (str): Path to the MOL2 file of the ligand (e.g., 'selected/1bai/1bai_ligand.mol2')
    - output_base_dir (str): Path to the base directory where results will be stored
    """

    # 1) Collect all XYZ files from chain_dir
    xyz_files = glob.glob(os.path.join(chain_dir, "*.xyz"))
    if not xyz_files:
        print(f"No .xyz files found in {chain_dir}")
        return

    # 2) Create the base output directory if it doesn't exist
    os.makedirs(output_base_dir, exist_ok=True)

    # 3) Loop over each XYZ file
    for xyz_file in xyz_files:
        # Extract a name for subdirectories / file naming
        xyz_name = os.path.splitext(os.path.basename(xyz_file))[0]

        # Create an output folder for this particular XYZ docking
        this_output_dir = os.path.join(output_base_dir, xyz_name)
        os.makedirs(this_output_dir, exist_ok=True)

        print(f"\n=== Processing {xyz_name} ===")

        # 4) Create a fresh Docking instance, specifying the output directory
        docking = Docking(output_dir=this_output_dir)

        # 5) Read the XYZ file
        docking.read_xyz(xyz_file)

        # 6) Adjust scale (optional: you can set a different target distance if desired)
        docking.adjust_scale(target_distance=3.8)

        # 7) Generate the full-atom PDB model
        #    Use a unique sequence name (e.g., xyz_name) to distinguish each model
        full_model_pdb = docking.generate_full_model(chain_id='A', sequence_name=xyz_name)
        if not full_model_pdb:
            print(f"Skipping docking for {xyz_file} because model generation failed.")
            continue

        # 8) Convert the full-atom PDB model to PDBQT for use as the receptor
        receptor_pdbqt = docking.convert_to_pdbqt(full_model_pdb)
        if not receptor_pdbqt:
            print(f"Skipping docking for {xyz_file} because receptor PDBQT conversion failed.")
            continue

        # 9) Run docking with AutoDock Vina
        #    - The ligand is the same for all runs: ligand_mol2
        #    - If you need a specific box size or exhaustiveness, adjust run_vina arguments
        docking_output_file, docking_log_file = docking.run_vina(
            receptor_pdbqt=receptor_pdbqt,
            ligand_mol2=ligand_mol2,
            seed=2,  # or another integer
            box_size=(18, 18, 18)
        )
        if docking_output_file and docking_log_file:
            print(f"Docking complete for {xyz_name}.")
            # 10) Optionally parse docking results (scores)
            scores = docking.parse_docking_results(docking_log_file)
            print(f"Docking scores for {xyz_name}: {scores}")
        else:
            print(f"Docking failed for {xyz_name}.")


if __name__ == "__main__":
    # Example usage:
    chain_5_directory = "grouped_result/chain_5/1fkf"  # e.g., "grouped_result/chain_5"
    ligand_file = "selected/1bai/1bai_ligand.mol2"  # e.g., "selected/1bai/1bai_ligand.mol2"
    results_directory = "docking_results/chain_5"  # e.g., "docking_results/chain_5"

    batch_docking_for_xyz(chain_5_directory, ligand_file, results_directory)







