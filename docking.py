# --*-- conding:utf-8 --*--
# @Time : 3/20/25 8:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking.py

from docking import Docking
import os

if __name__ == "__main__":

    xyz_path = "1fkf.xyz"  # 请修改为实际路径
    docking_folder = "./pdbqt_protein/docking_output"  # 指定保存 PDBQT 的文件夹
    recon = Docking(xyz_path, docking_folder=docking_folder)
    recon.run_pipeline()







