# --*-- conding:utf-8 --*--
# @Time : 3/21/25 12:57 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Docking.py


import os
import numpy as np
import shutil

# Modeller相关
from modeller import environ, log
from modeller.automodel import automodel, assess
from modeller import alignment, model

class Docking:
    """
    使用MODELLER从仅含Cα的xyz文件构建满原子PDB的示例类。

    与原先类保持方法和流程一致:
      1) read_xyz()         # 读取xyz(只含Cα)
      2) adjust_scale()     # 缩放Cα
      3) write_ca_pdb()     # 写出仅含Cα的pdb
      4) prepare_alignment()# 写出alignment文件(模板 + 目标序列)
      5) generate_full_model() # 执行MODELLER的automodel建模

    最终可以把得到的full_model重命名或复制为 self.full_model_filename。
    """

    def __init__(self, xyz_file):
        """
        参数
        ----
        xyz_file : str
            仅含Cα坐标的xyz文件路径
        """
        self.xyz_file = xyz_file
        self.sequence = []
        self.coordinates = []
        self.scaled_coordinates = []

        # 确定 xyz_file 所在目录，后面写文件都放在同一目录下
        self.base_dir = os.path.dirname(os.path.abspath(self.xyz_file))
        print("XYZ base directory:", self.base_dir)

        # 基于 xyz_file 命名一些中间文件
        self.ca_pdb_filename = os.path.join(self.base_dir, 'ca_model.pdb')  # 只含Cα的 PDB
        self.alignment_filename = os.path.join(self.base_dir, 'protein.ali')  # alignment 文件
        self.full_model_filename = os.path.join(self.base_dir, 'full_model.pdb')  # 期望得到的完整PDB

        # 这里根据 xyz_file 的倒数第5个字符来做 chain_id，原先你是这么写的
        # 若不符合你的命名习惯，可自行修改
        self.chain_id = self.xyz_file[-5]

    def read_xyz(self):
        """
        从 xyz 文件中读取残基的一字母名称与Cα坐标。
        注意：本方法跳过前两行（通常第一行是原子数，第二行是注释）。
        """
        with open(self.xyz_file, 'r') as f:
            lines = f.readlines()[2:]  # 跳过前两行（有时第一行是原子数，第二行是空行/注释）
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
        计算当前Cα的平均相邻距离，并用target_distance(默认3.8Å)做整体缩放。
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
        self.scaled_coordinates = [
            (x * scale_factor, y * scale_factor, z * scale_factor)
            for x, y, z in self.coordinates
        ]
        print(f"Scale factor applied: {scale_factor:.4f}")

    def write_ca_pdb(self):
        """
        将缩放后的Cα坐标写成一个仅含Cα的PDB文件(ATOM行)。
        """
        # 一字母 -> 三字母映射
        res_name_map = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
            'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }

        with open(self.ca_pdb_filename, 'w') as f:
            for i, (res_name, (x, y, z)) in enumerate(zip(self.sequence, self.scaled_coordinates), start=1):
                # 三字母
                res_name_3 = res_name_map.get(res_name.upper(), 'UNK')
                # 这里的 chain_id = self.chain_id, 残基编号从1开始
                # 原子名称我们统一写 CA（仅Cα）
                f.write(
                    f"ATOM  {i:5d}  CA  {res_name_3} {self.chain_id}{i:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n"
                )
            f.write("END\n")
        print(f"PDB file '{self.ca_pdb_filename}' has been written.")

    def prepare_alignment(self, sequence_name='protein_full'):
        """
        写出 alignment 文件，以用于 MODELLER 的 automodel。

        思路：把仅含Cα的ca_model当模板；把目标序列命名为 protein_full。
        让这两者在对齐文件中使用相同的氨基酸序列 (只是一字母拼起来)，加个'*'结束。
        """
        sequence_str = ''.join(self.sequence) + '*'
        start_residue = 1
        end_residue = len(self.sequence)

        with open(self.alignment_filename, 'w') as f:
            # 模板（来自ca_model）
            f.write(">P1;ca_model\n")
            # structureX:ca_model:<start>:<chain>:<end>:<chain>::::
            # 这里 chain 都用 self.chain_id，res编号从start_residue到end_residue
            f.write(f"structureX:ca_model:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n\n")

            # 目标（full model）
            f.write(f">P1;{sequence_name}\n")
            f.write(f"sequence:{sequence_name}:{start_residue}:{self.chain_id}:{end_residue}:{self.chain_id}::::\n")
            f.write(sequence_str + "\n")

        print(f"Alignment file '{self.alignment_filename}' has been prepared.")

    def generate_full_model(self):
        """
        用 MODELLER 的 automodel 从仅含Cα的ca_model来构建完整蛋白。
        将结果输出到默认命名 (protein_full.B9999000x.pdb)，
        然后复制/重命名为 self.full_model_filename。
        """

        # 关闭或开启 Modeller 日志
        log.none()
        env = environ()

        # (可选) 如果需要载入自带的top_heav.lib等，可以：
        # env.libs.topology.read(file='$(LIB)/top_heav.lib')
        # env.libs.parameters.read(file='$(LIB)/par.lib')
        #
        # 同时如果需要让Modeller能找到ca_model.pdb, alignment等文件：
        env.io.atom_files_directory = [self.base_dir]

        # 先做对齐
        aln = alignment(env)
        # 读入仅含Cα的 pdb，作为 knowns='ca_model'
        mdl = model(env, file=self.ca_pdb_filename,
                    model_segment=(f'FIRST:{self.chain_id}', f'LAST:{self.chain_id}'))
        aln.append_model(mdl, align_codes='ca_model', atom_files=self.ca_pdb_filename)

        # 将 alignment 文件读入(包含 'protein_full' 目标序列)
        aln.append(file=self.alignment_filename, align_codes='protein_full')
        # 两个序列做一次二级对齐
        aln.align2d()

        # 定义 automodel 类
        class MyModel(automodel):
            def special_patches(self, aln):
                # 如果需要修饰链ID、指定环区补建等可在此处理
                self.rename_segments(segment_ids=['A'])

        # 构建
        a = MyModel(
            env,
            alnfile=self.alignment_filename,
            knowns='ca_model',
            sequence='protein_full',
            assess_methods=(assess.DOPE, assess.GA341),
        )
        a.starting_model = 1
        a.ending_model = 1  # 只做1个模型
        a.make()

        # a.outputs 里存储了建模结果信息，比如最佳模型名字
        # 形如 [{'name': 'protein_full.B99990001.pdb', 'mdl': <modeller.model.Model object at 0x...>}]
        best_model_path = a.outputs[0]['name']
        print("MODELLER built a model:", best_model_path)

        # 将生成的PDB复制(或重命名)为 self.full_model_filename
        shutil.copyfile(best_model_path, self.full_model_filename)
        print(f"Final full model has been saved as '{self.full_model_filename}'.")

    def run_pipeline(self):
        """
        可选：把所有步骤封装到一个一键执行的方法，用户只需调用 run_pipeline() 即可生成full model。
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
        print("==== DONE ====")


















