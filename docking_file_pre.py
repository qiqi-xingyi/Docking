# --*-- conding:utf-8 --*--
# @Time : 3/20/25 8:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking_file_pre.py

from docking import Fileprepare
import os

# 示例使用：
if __name__ == "__main__":
    # 处理蛋白：从 XYZ 到 PDBQT
    xyz_path = "1fkf.xyz"  # 请修改为实际路径
    docking_folder = "./pdbqt"  # 指定保存最终 PDBQT 的文件夹
    docking_obj = Fileprepare(xyz_path, docking_folder=docking_folder)
    docking_obj.run_pipeline()

    # 处理配体：从 MOL2 到 pdbqt
    ligand_mol2_path = "selected/1fkf/1fkf_ligand.mol2"  # 请修改为实际路径
    # 此处可以复用 docking_obj 的 docking_folder 来保存配体的 pdbqt 文件
    ligand_pdbqt = docking_obj.prepare_ligand_pdbqt(ligand_mol2_path, output_folder=docking_folder)
    print("Final ligand pdbqt file:", ligand_pdbqt)







