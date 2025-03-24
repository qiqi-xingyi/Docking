# --*-- conding:utf-8 --*--
# @Time : 3/20/25 8:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking_file_pre.py

import os

from docking import Fileprepare  # 假设Docking类和Mol2Translator类都定义在Docking.py里

# 1) 设置主要路径
grouped_result_dir = "./grouped_result"  # 包含1fkf、3ckz等子目录，每个目录有多个xyz文件
selected_dir = "./selected"  # 包含1fkf、3ckz等子目录，每个目录有ligand.mol2
output_protein_dir = "./proteins_pdbqt"  # 蛋白PDBQT输出根目录
output_ligand_dir = "./ligands_pdbqt"  # 配体PDBQT输出根目录

# 2) 确保输出目录存在
os.makedirs(output_protein_dir, exist_ok=True)
os.makedirs(output_ligand_dir, exist_ok=True)


def extract_protein_id(xyz_filename):
    """
    从xyz文件名中提取蛋白ID的示例函数。
    假设文件名可能是 '1fkf.xyz' 或 '1fkf_top_1.xyz' 等，
    我们先去掉'.xyz'，再按'_'分割，取第一个部分作为蛋白ID。
    """
    base = os.path.splitext(xyz_filename)[0]  # e.g. "1fkf_top_1" 或 "1fkf"
    parts = base.split('_')  # e.g. ["1fkf", "top", "1"] 或 ["1fkf"]
    return parts[0]  # -> "1fkf"

if __name__ == '__main__':

    # 3) 递归或逐层遍历 grouped_result_dir 下的所有子目录
    for protein_dir in os.listdir(grouped_result_dir):
        protein_path = os.path.join(grouped_result_dir, protein_dir)
        # 跳过非目录
        if not os.path.isdir(protein_path):
            continue

        # 遍历当前目录下的 xyz 文件
        for file_name in os.listdir(protein_path):
            if not file_name.endswith(".xyz"):
                continue

            # ========== 跳过与 _top_1.xyz 重复的“基本文件” ==========
            # 如果文件名不包含 "_top_" (比如 "1fkf.xyz"),
            # 且同目录下存在 "1fkf_top_1.xyz"，就跳过这个文件。
            if "_top_" not in file_name:
                # 提取ID，比如 "1fkf.xyz" -> "1fkf"
                protein_id_temp = extract_protein_id(file_name)  # "1fkf"
                # 构造对应的 "_top_1" 文件名： "1fkf_top_1.xyz"
                top_1_filename = f"{protein_id_temp}_top_1.xyz"
                top_1_path = os.path.join(protein_path, top_1_filename)
                if os.path.exists(top_1_path):
                    print(f"Skipping {file_name} because {top_1_filename} exists (they are duplicates).")
                    continue

            # ========== 处理蛋白 (XYZ -> PDBQT) ==========
            xyz_path = os.path.join(protein_path, file_name)
            protein_id = extract_protein_id(file_name)

            # 创建输出目录
            protein_output_dir = os.path.join(output_protein_dir, protein_id)
            os.makedirs(protein_output_dir, exist_ok=True)

            print(f"\n=== Processing protein '{protein_id}' from '{file_name}' ===")
            docking_obj = Fileprepare(
                xyz_file=xyz_path,
                docking_folder=protein_output_dir
            )
            docking_obj.run_pipeline()

            # ========== 处理配体 (MOL2 -> PDBQT) ==========
            # 假设配体文件在 selected_dir/<protein_id>/<protein_id>_ligand.mol2
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









