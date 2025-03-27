# --*-- conding:utf-8 --*--
# @Time : 3/25/25 6:47 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : cif_prepare.py

import os
from docking import Cifprepare

def process_all_cif_files(source_dir, target_dir):

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)


    for folder_name in os.listdir(source_dir):
        subdir_path = os.path.join(source_dir, folder_name)
        if os.path.isdir(subdir_path):

            last_cif_file = None
            for fname in os.listdir(subdir_path):
                if fname.endswith("4.cif"):
                    last_cif_file = fname

                    break

            if last_cif_file:

                cif_path = os.path.join(subdir_path, last_cif_file)

                dest_subdir = os.path.join(target_dir, folder_name)
                os.makedirs(dest_subdir, exist_ok=True)

                print(f"Start {cif_path} ...")
                processor = Cifprepare(cif_file=cif_path, output_folder=dest_subdir)

                processor.run_pipeline(translate=True)

                print(f"Completed：{cif_path} => {dest_subdir}")
            else:
                print(f"Warning： {folder_name} don't have file end with '4.cif', pass.")

if __name__ == "__main__":

    source_folder = "af3/result_af3_2"
    target_folder = "af3_pdbqt"

    process_all_cif_files(source_folder, target_folder)
