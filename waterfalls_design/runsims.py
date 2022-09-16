import csv
import numpy as np
import os
import math
import sys
from os import listdir
import time
import shutil

files_to_run = ['waterfalls_20','waterfalls_21']


for f in files_to_run:
    gen_path = os.getcwd() + "/" + f + "/mesh_generator.py"
    os.system('salome shell ' + gen_path)
    shutil.move(os.getcwd() + "/mesh.unv", os.getcwd() + "/" + f + "/mesh.unv")

for f in files_to_run:
    gen_path = os.getcwd() + "/" + f + "/mesh.unv"
    os.system('ElmerGrid 8 2 ' + gen_path)

for f in files_to_run:
    gen_path = os.getcwd() + "/" + f + "/sif_generator_from_code.py"
    os.system('python3 ' + gen_path)

for f in files_to_run:
    gen_path = os.getcwd() + "/" + f + "/mesh/mesh.sif"
    os.system('ElmerSolver ' + gen_path)
