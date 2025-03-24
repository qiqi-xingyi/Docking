# --*-- conding:utf-8 --*--
# @Time : 3/20/25 8:43 PM
# @Author : Yuqi Zhang
# @Email : yzhan135@kent.edu
# @File : docking.py

from docking import Docking
import os

if __name__ == "__main__":

    xyz_path = "3ckz.xyz"
    recon = Docking(xyz_path)
    recon.run_pipeline()







