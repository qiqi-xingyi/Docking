# --*-- conding:utf-8 --*--
# @Time : 3/25/25 6:47 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : cif_prepare.py


from docking import Cifprepare


if __name__ == "__main__":

    cif_file_path = "example.cif"
    docking_dir = "docking_results"


    processor = Cifprepare(cif_file_path, output_folder=docking_dir)
    processor.run_pipeline(translate=True)
