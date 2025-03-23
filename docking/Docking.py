# --*-- conding:utf-8 --*--
# @Time : 3/21/25 12:57 AM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : Docking.py

import os
import subprocess
import tempfile
import numpy as np
from modeller import environ, automodel, alignment, model, log
from modeller.automodel import assess


class Docking:
    def __init__(self, output_dir=None):
        """
        Initialize the Docking pipeline.
        If output_dir is not provided, the system temporary directory is used.
        """
        if output_dir is None:
            self.output_dir = tempfile.gettempdir()
        else:
            self.output_dir = output_dir
            os.makedirs(self.output_dir, exist_ok=True)

        # Modeling attributes
        self.sequence = []
        self.coordinates = []
        self.scaled_coordinates = []
        self.xyz_file = None

    # ------------------ Modeling Methods ------------------ #
    def read_xyz(self, xyz_file):
        """
        Read an XYZ file (ignoring the first two lines) and store the sequence and coordinates.
        """
        self.xyz_file = xyz_file
        with open(xyz_file, 'r') as f:
            lines = f.readlines()[2:]
        for line in lines:
            parts = line.strip().split()
            if len(parts) == 4:
                res_name = parts[0]
                x, y, z = map(float, parts[1:])
                self.sequence.append(res_name)
                self.coordinates.append((x, y, z))
            else:
                print(f"Warning: Malformed line skipped: {line.strip()}")
        return self.sequence, self.coordinates

    def adjust_scale(self, target_distance=3.8):
        """
        Scale coordinates so that the average distance between consecutive residues is target_distance.
        """

        def calculate_average_distance(coords):
            distances = []
            for i in range(len(coords) - 1):
                x1, y1, z1 = coords[i]
                x2, y2, z2 = coords[i + 1]
                dist = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
                distances.append(dist)
            return np.mean(distances)

        current_avg = calculate_average_distance(self.coordinates)
        scale_factor = target_distance / current_avg
        self.scaled_coordinates = [(x * scale_factor, y * scale_factor, z * scale_factor)
                                   for x, y, z in self.coordinates]
        print(f"Scale factor applied: {scale_factor:.4f}")
        return scale_factor

    def generate_ca_pdb_content(self, chain_id='A'):
        """
        Generate the CA model content in PDB format as a string.
        """
        res_name_map = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
            'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
            'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
            'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
            'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }
        pdb_lines = []
        for i, (res, (x, y, z)) in enumerate(zip(self.sequence, self.scaled_coordinates), start=1):
            res_3 = res_name_map.get(res.upper(), 'UNK')
            pdb_lines.append(
                f"ATOM  {i:5d}  CA  {res_3} {chain_id}{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
            )
        pdb_lines.append("END")
        return "\n".join(pdb_lines)

    def generate_alignment_content(self, chain_id='A', sequence_name='protein_full'):
        """
        Generate the alignment file content (as a string) for Modeller.
        """
        sequence_str = ''.join(self.sequence) + '*'
        start_res = 1
        end_res = len(self.sequence)
        lines = []
        # Template entry (CA model)
        lines.append(">P1;ca_model")
        lines.append(f"structureX:ca_model:{start_res}:{chain_id}:{end_res}:{chain_id}::::")
        lines.append(sequence_str)
        lines.append("")
        # Target entry (full model)
        lines.append(f">P1;{sequence_name}")
        lines.append(f"sequence:{sequence_name}:{start_res}:{chain_id}:{end_res}:{chain_id}::::")
        lines.append(sequence_str)
        return "\n".join(lines)

    def generate_full_model(self, chain_id='A', sequence_name='protein_full'):
        """
        Generate a full atomic model using Modeller.
        Intermediate CA PDB and alignment files are saved as temporary files and deleted afterward.
        """
        ca_pdb_content = self.generate_ca_pdb_content(chain_id=chain_id)
        aln_content = self.generate_alignment_content(chain_id=chain_id, sequence_name=sequence_name)

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.pdb', delete=False, dir=self.output_dir) as ca_file:
            ca_file.write(ca_pdb_content)
            ca_file_path = ca_file.name

        with tempfile.NamedTemporaryFile(mode='w+', suffix='.ali', delete=False, dir=self.output_dir) as aln_file:
            aln_file.write(aln_content)
            aln_file_path = aln_file.name


        log.none()
        env = environ()
        env.io.atom_files_directory = [self.output_dir]
        aln = alignment(env)
        mdl = model(env, file=ca_file_path, model_segment=(f'FIRST:{chain_id}', f'LAST:{chain_id}'))
        aln.append_model(mdl, align_codes='ca_model', atom_files=ca_file_path)
        aln.append(file=aln_file_path, align_codes=sequence_name)
        aln.align2d()

        class MyModel(automodel):
            def special_patches(self, aln):
                self.rename_segments(segment_ids=[chain_id])

        a = MyModel(env, alnfile=aln_file_path, knowns='ca_model', sequence=sequence_name,
                    assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model = 1
        a.make()

        full_model_path = os.path.join(self.output_dir, f"{sequence_name}.B99990001.pdb")
        os.remove(ca_file_path)
        os.remove(aln_file_path)

        print(f"Full model generated: {full_model_path}")
        return full_model_path

    # ------------------ Docking Methods ------------------ #
    def _is_tool_available(self, tool_name):
        """
        Check if a tool (e.g., obabel or vina) is available on the system.
        """
        from shutil import which
        return which(tool_name) is not None

    def convert_mol2_to_pdbqt(self, mol2_file, pdbqt_file):
        """
        Convert a MOL2 file to PDBQT format using Open Babel.
        """
        if not self._is_tool_available("obabel"):
            raise EnvironmentError("Open Babel is not installed or not in PATH.")
        command = ["obabel", mol2_file, "-O", pdbqt_file, "--partialcharge", "gasteiger", "-h"]
        subprocess.run(command, check=True)
        print(f"Converted {mol2_file} to {pdbqt_file}.")

    def calculate_ligand_center_of_mass(self, mol2_file):
        """
        Calculate the center of mass of the ligand from a MOL2 file.
        """
        atom_coords = []
        with open(mol2_file, 'r') as file:
            atom_section = False
            for line in file:
                if line.startswith('@<TRIPOS>ATOM'):
                    atom_section = True
                    continue
                elif line.startswith('@<TRIPOS>') and atom_section:
                    break
                if atom_section:
                    parts = line.split()
                    if len(parts) >= 5:
                        x, y, z = map(float, parts[2:5])
                        atom_coords.append((x, y, z))
        if atom_coords:
            x_coords, y_coords, z_coords = zip(*atom_coords)
            center_x = sum(x_coords) / len(x_coords)
            center_y = sum(y_coords) / len(y_coords)
            center_z = sum(z_coords) / len(z_coords)
            print(f"Calculated ligand center of mass: X={center_x:.3f}, Y={center_y:.3f}, Z={center_z:.3f}")
            return center_x, center_y, center_z
        else:
            raise ValueError("No atom coordinates found in the MOL2 file.")

    def run_vina(self, receptor_pdbqt, ligand_mol2, seed=2, box_size=(18, 18, 18)):
        """
        Run AutoDock Vina docking.
        Converts the ligand from MOL2 to PDBQT, calculates the center of mass, and then runs docking.
        """
        # Convert ligand from MOL2 to PDBQT
        ligand_pdbqt = os.path.join(self.output_dir, "ligand.pdbqt")
        self.convert_mol2_to_pdbqt(ligand_mol2, ligand_pdbqt)

        # Calculate center of mass for the docking box from the ligand MOL2 file
        center_x, center_y, center_z = self.calculate_ligand_center_of_mass(ligand_mol2)

        # Define output files
        output_file = os.path.join(self.output_dir, "docking_output.pdbqt")
        log_file = os.path.join(self.output_dir, "docking_log.txt")

        size_x, size_y, size_z = box_size
        command = [
            "vina",
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--out", output_file,
            "--center_x", str(center_x), "--center_y", str(center_y), "--center_z", str(center_z),
            "--size_x", str(size_x), "--size_y", str(size_y), "--size_z", str(size_z),
            "--exhaustiveness", "16",
            "--seed", str(seed)
        ]
        print("Running AutoDock Vina docking...")
        try:
            with open(log_file, "w") as log:
                subprocess.run(command, stdout=log, stderr=log, check=True)
        except subprocess.CalledProcessError as e:
            print("Error: AutoDock Vina docking failed.")
            print("Command:", " ".join(command))
            print("Exit status:", e.returncode)
            return None, None
        print(f"Docking complete. Results saved in {output_file}")
        return output_file, log_file

    def parse_docking_results(self, log_file):
        """
        Parse the AutoDock Vina log file to extract docking scores.
        """
        scores = []
        with open(log_file, "r") as file:
            for line in file:
                if "REMARK VINA RESULT:" in line:
                    parts = line.strip().split()
                    # Binding affinity is typically the fourth element
                    score = float(parts[3])
                    scores.append(score)
        print(f"Extracted docking scores: {scores}")
        return scores

    # ------------------ MOL2 Translation Methods ------------------ #
    def translate_mol2(self, input_mol2, output_mol2):
        """
        Translate the molecule in a MOL2 file so that its geometric center is at the origin.
        """
        atoms = []
        # Read the MOL2 file and extract atom information
        with open(input_mol2, 'r') as file:
            lines = file.readlines()

        atom_section = False
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                atom_section = True
                continue
            elif line.startswith("@<TRIPOS>") and atom_section:
                break
            if atom_section:
                parts = line.split()
                if len(parts) < 6:
                    continue
                atom_id = parts[0]
                atom_name = parts[1]
                try:
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                except ValueError:
                    raise ValueError(f"Invalid coordinates for atom {atom_id} in {input_mol2}.")
                atom_type = parts[5]
                additional = parts[6:] if len(parts) > 6 else []
                atoms.append({
                    'atom_id': atom_id,
                    'atom_name': atom_name,
                    'x': x,
                    'y': y,
                    'z': z,
                    'atom_type': atom_type,
                    'additional': additional
                })

        if not atoms:
            raise ValueError(f"No atoms found in {input_mol2}.")

        # Calculate geometric center
        num_atoms = len(atoms)
        center_x = sum(atom['x'] for atom in atoms) / num_atoms
        center_y = sum(atom['y'] for atom in atoms) / num_atoms
        center_z = sum(atom['z'] for atom in atoms) / num_atoms
        print(f"Geometric center: X={center_x:.3f}, Y={center_y:.3f}, Z={center_z:.3f}")

        # Translate atoms so that the center is at the origin
        for atom in atoms:
            atom['x'] -= center_x
            atom['y'] -= center_y
            atom['z'] -= center_z

        # Write the translated MOL2 file while preserving original structure
        with open(output_mol2, 'w') as outfile:
            atom_section_flag = False
            for line in lines:
                stripped = line.strip()
                if stripped.startswith("@<TRIPOS>ATOM"):
                    atom_section_flag = True
                    outfile.write(line)
                    continue
                elif stripped.startswith("@<TRIPOS>") and atom_section_flag:
                    atom_section_flag = False
                    outfile.write(line)
                    continue

                if atom_section_flag:
                    parts = line.split()
                    if len(parts) < 6:
                        outfile.write(line)
                        continue
                    atom_id = parts[0]
                    atom_data = next((a for a in atoms if a['atom_id'] == atom_id), None)
                    if atom_data is None:
                        outfile.write(line)
                        continue
                    new_line = f"{atom_data['atom_id']:>6} {atom_data['atom_name']:<10} {atom_data['x']:>8.3f} {atom_data['y']:>8.3f} {atom_data['z']:>8.3f} {atom_data['atom_type']}"
                    if atom_data['additional']:
                        new_line += " " + " ".join(atom_data['additional'])
                    new_line += "\n"
                    outfile.write(new_line)
                else:
                    outfile.write(line)
        print(f"Translated MOL2 file saved as {output_mol2}.")

    def convert_to_pdbqt(self, input_file, output_file=None):
        """
        Convert a PDB or MOL2 file to PDBQT format using Open Babel.
        Automatically detects the file type from extension (or let Open Babel handle it).
        If output_file is not provided, create one in self.output_dir.
        """
        if not self._is_tool_available("obabel"):
            raise EnvironmentError("Open Babel is not installed or not in PATH.")

        if output_file is None:
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            output_file = os.path.join(self.output_dir, f"{base_name}.pdbqt")

        command = [
            "obabel", input_file,
            "-O", output_file,
            "--partialcharge", "gasteiger",
            "-h"  # add hydrogens
        ]
        subprocess.run(command, check=True)
        print(f"Converted {input_file} to {output_file}.")
        return output_file



