# --*-- conding:utf-8 --*--
# @Time : 3/21/25 12:57 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Docking.py


import os
import numpy as np
import subprocess
import shutil
from modeller import environ, log
from modeller.automodel import automodel, assess
from modeller import alignment, model
from Bio.PDB import PDBParser, PDBIO, MMCIFParser


class Docking:
    """
    使用MODELLER从仅含 Cα 的 XYZ 文件构建满原子 PDB，并进一步添加氢、配平电荷，
    转换为 PDBQT 文件。最终的 PDBQT 文件将保存至指定文件夹中。

    流程：
      1) read_xyz()           : 读取 XYZ 文件（仅含 Cα 坐标及残基一字母代码）
      2) adjust_scale()       : 缩放坐标，使相邻 Cα 距离为 3.8 Å
      3) write_ca_pdb()       : 写出仅含 Cα 的 PDB 文件（作为 MODELLER 模板）
      4) prepare_alignment()  : 生成用于建模的 alignment 文件
      5) generate_full_model(): 使用 automodel 构建满原子模型（完整 PDB）
      6) prepare_pdbqt()      : 添加氢原子、配平电荷并转换为 PDBQT 文件
      7) run_pipeline()       : 一键运行上述全部步骤
    """

    def __init__(self, xyz_file, docking_folder=None):
        """
        Parameters:
        -----------
        xyz_file : str
            包含仅 Cα 坐标和残基一字母代码的 XYZ 文件路径
        docking_folder : str, optional
            指定保存最终 PDBQT 文件的文件夹，默认为 xyz 文件所在目录
        """
        self.xyz_file = xyz_file
        self.sequence = []
        self.coordinates = []
        self.scaled_coordinates = []

        # 确定 xyz 文件所在目录
        self.base_dir = os.path.dirname(os.path.abspath(self.xyz_file))
        print("XYZ base directory:", self.base_dir)

        # 根据 xyz 文件生成中间文件名
        self.ca_pdb_filename = os.path.join(self.base_dir, 'ca_model.pdb')  # 仅含 Cα 的 PDB
        self.alignment_filename = os.path.join(self.base_dir, 'protein.ali')  # alignment 文件
        self.full_model_filename = os.path.join(self.base_dir, 'full_model.pdb')  # MODELLER生成的满原子PDB

        # 最终的 PDBQT 文件名（初步在 base_dir 中生成，后续可移动到 docking_folder）
        self.pdbqt_filename = self.full_model_filename.replace(".pdb", ".pdbqt")

        # 如果提供了 docking_folder，则保存最终 PDBQT 到该目录，否则默认为 base_dir
        if docking_folder is None:
            self.docking_folder = self.base_dir
        else:
            self.docking_folder = docking_folder
            os.makedirs(self.docking_folder, exist_ok=True)

        # 这里根据 xyz 文件名的倒数第5个字符设定 chain_id，若不适用请自行修改
        self.chain_id = self.xyz_file[-5]

    def read_xyz(self):
        """读取 XYZ 文件，提取残基一字母代码和 Cα 坐标。跳过前两行。"""
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()[2:]  # 跳过前两行（通常为原子数与注释）
            for line in lines:
                parts = line.strip().split()
                if len(parts) == 4:
                    res_name = parts[0]
                    x, y, z = map(float, parts[1:])
                    self.sequence.append(res_name)
                    self.coordinates.append((x, y, z))
                else:
                    print(f"Warning: Line '{line.strip()}' is malformed and will be skipped.")

    def adjust_scale(self, target_distance=3.8):
        """
        根据相邻 Cα 之间的平均距离计算缩放因子，使其调整为 target_distance（默认 3.8 Å）。
        """

        def calculate_average_distance(coords):
            distances = []
            for i in range(len(coords) - 1):
                x1, y1, z1 = coords[i]
                x2, y2, z2 = coords[i + 1]
                distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
                distances.append(distance)
            return np.mean(distances)

        if len(self.coordinates) < 2:
            print("Not enough residues to compute distances.")
            self.scaled_coordinates = self.coordinates[:]
            return

        current_avg_distance = calculate_average_distance(self.coordinates)
        scale_factor = target_distance / current_avg_distance
        self.scaled_coordinates = [(x * scale_factor, y * scale_factor, z * scale_factor)
                                   for x, y, z in self.coordinates]
        print(f"Scale factor applied: {scale_factor:.4f}")

    def write_ca_pdb(self):
        """
        将缩放后的 Cα 坐标写入 PDB 文件，仅包含 CA 原子信息，供 MODELLER 补全用。
        """
        res_name_map = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
            'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }
        with open(self.ca_pdb_filename, 'w') as f:
            for i, (res_name, (x, y, z)) in enumerate(zip(self.sequence, self.scaled_coordinates), start=1):
                res_name_3 = res_name_map.get(res_name.upper(), 'UNK')
                f.write(f"ATOM  {i:5d}  CA  {res_name_3} {self.chain_id}{i:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
            f.write("END\n")
        print(f"PDB file '{self.ca_pdb_filename}' has been written.")

    def prepare_alignment(self, sequence_name='protein_full'):
        """
        根据残基序列准备 alignment 文件，其中模板为 ca_model，目标为 protein_full。
        注意序列末尾添加 '*' 以结束序列。
        """
        sequence_str = ''.join(self.sequence) + '*'
        start_residue = 1
        end_residue = len(self.sequence)

        with open(self.alignment_filename, 'w') as f:
            # 模板部分（来自 ca_model）
            f.write(">P1;ca_model\n")
            f.write(f"structureX:ca_model:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n\n")
            # 目标部分（full model）
            f.write(f">P1;{sequence_name}\n")
            f.write(f"sequence:{sequence_name}:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n")
        print(f"Alignment file '{self.alignment_filename}' has been prepared.")

    def generate_full_model(self):
        """
        使用 MODELLER 的 automodel 方法从仅含 Cα 的 ca_model 补全构建满原子结构，
        并将生成的 PDB 复制为 full_model.pdb。
        """
        log.none()  # 如需详细输出，可使用 log.verbose()
        env = environ()
        env.io.atom_files_directory = [self.base_dir]

        # 读入仅含 Cα 的 PDB 作为模板
        aln = alignment(env)
        mdl = model(env, file=self.ca_pdb_filename,
                    model_segment=(f'FIRST:{self.chain_id}', f'LAST:{self.chain_id}'))
        aln.append_model(mdl, align_codes='ca_model', atom_files=self.ca_pdb_filename)
        aln.append(file=self.alignment_filename, align_codes='protein_full')
        aln.align2d()

        # 定义 automodel 子类，可在 special_patches 中修改链ID或其它修饰操作
        class MyModel(automodel):
            def special_patches(self, aln):
                self.rename_segments(segment_ids=['A'])

        a = MyModel(env,
                    alnfile=self.alignment_filename,
                    knowns='ca_model',
                    sequence='protein_full',
                    assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        best_model_path = a.outputs[0]['name']
        print("MODELLER built a model:", best_model_path)
        shutil.copyfile(best_model_path, self.full_model_filename)
        print(f"Final full model has been saved as '{self.full_model_filename}'.")

    # 以下为添加氢、配平电荷、转换成 PDBQT 的新功能
    def _is_tool_available(self, tool_name):
        """检查系统中是否有指定工具（如 obabel）。"""
        return subprocess.call(f"type {tool_name}", shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

    def calculate_center_of_mass(self, structure_file):
        """
        计算指定 PDB 文件中所有原子的中心质心。
        """
        # 这里仅支持 PDB 格式（若为 CIF，可引入 MMCIFParser）
        parser = PDBParser(QUIET=True)
        try:
            structure = parser.get_structure('structure', structure_file)
        except Exception as e:
            raise ValueError(f"Failed to parse file {structure_file}: {e}")
        atom_coords = []
        for atom in structure.get_atoms():
            atom_coords.append(atom.get_coord())
        if not atom_coords:
            raise ValueError(f"No atom coordinates found in file {structure_file}.")
        atom_coords = np.array(atom_coords)
        center_of_mass = atom_coords.mean(axis=0)
        print(f"Center of mass: X={center_of_mass[0]:.3f}, Y={center_of_mass[1]:.3f}, Z={center_of_mass[2]:.3f}")
        return center_of_mass

    def translate_to_origin(self, input_file, output_file):
        """
        将分子平移，使其中心质心移动到原点，然后保存为新文件。
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', input_file)
        center_of_mass = self.calculate_center_of_mass(input_file)
        translation = -center_of_mass  # 计算平移向量
        for atom in structure.get_atoms():
            atom.set_coord(atom.get_coord() + translation)
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)
        print(f"Molecule has been translated to the origin and saved as {output_file}")

    def _convert_to_pdbqt_with_obabel(self, input_file, output_pdbqt):
        """使用 Open Babel 将 PDB 添加氢、配平电荷并转换为 PDBQT格式。"""
        command = [
            "obabel", input_file, "-O", output_pdbqt,
            "--partialcharge", "gasteiger", "-h", "-xr"
        ]
        subprocess.run(command, check=True)

    def prepare_pdbqt(self, translate=True, output_translated_file=None):
        """
        准备 PDBQT 文件：
          - 检查 Open Babel 是否可用；
          - （可选）将结构平移至原点；
          - 调用 Open Babel 添加氢、计算 Gasteiger 电荷，并生成 PDBQT 文件。
        """
        if not self._is_tool_available("obabel"):
            raise EnvironmentError("Open Babel is not installed or not in PATH.")

        if translate:
            if output_translated_file is None:
                output_translated_file = self.full_model_filename.replace(".pdb", "_translated.pdb")
            print("Translating molecule to the origin...")
            self.translate_to_origin(self.full_model_filename, output_translated_file)
            input_for_conversion = output_translated_file
        else:
            input_for_conversion = self.full_model_filename

        print("Converting to PDBQT using Open Babel...")
        self._convert_to_pdbqt_with_obabel(input_for_conversion, self.pdbqt_filename)
        print(f"Conversion complete. Intermediate PDBQT file: {self.pdbqt_filename}")

        # 如果指定了不同的 docking_folder，则将 pdbqt 移动过去
        if os.path.abspath(self.docking_folder) != os.path.abspath(self.base_dir):
            final_pdbqt_path = os.path.join(self.docking_folder, os.path.basename(self.pdbqt_filename))
            shutil.move(self.pdbqt_filename, final_pdbqt_path)
            self.pdbqt_filename = final_pdbqt_path
            print(f"PDBQT file has been moved to {self.pdbqt_filename}")

    def run_pipeline(self):
        """
        一键运行从 XYZ 到 PDBQT 的全部流程。
        """
        print("==== 1) Reading XYZ ====")
        self.read_xyz()
        print("==== 2) Adjusting scale ====")
        self.adjust_scale()
        print("==== 3) Writing CA PDB ====")
        self.write_ca_pdb()
        print("==== 4) Preparing alignment ====")
        self.prepare_alignment(sequence_name='protein_full')
        print("==== 5) Generating full model ====")
        self.generate_full_model()
        print("==== 6) Preparing PDBQT file ====")
        self.prepare_pdbqt(translate=True)
        print("==== Pipeline DONE ====")
        print(f"Final docking file: {self.pdbqt_filename}")






















