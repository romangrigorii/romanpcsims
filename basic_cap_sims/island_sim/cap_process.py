#!/usr/bin/env python
# This file generates a series of FEM meshes which will then be used to compute capacitances for

import os
from os import listdir

indir = listdir("UNV/")
dirs = []
for f in indir:
    if not os.path.isfile("UNV/" + f):
        dirs.append(f)

fw = open("capacitor_vals_i.txt","w")

for f in dirs:

    fr = open("UNV/" + f + "/" + f + "_cap","r")
    line = fr.readline()
    entries = line.split()
    fw.write(entries[1] + '\n')
