#!/usr/bin/env python
# This file generates a series of FEM meshes which will then be used to compute capacitances for

import os
from os import listdir

indir = listdir("UNV/")
dirs = []
for f in indir:
    if not os.path.isfile("UNV/" + f):
        dirs.append(f)

for f in dirs:
    fw = open("UNV/" + f + "/mutual_cap_vals.txt","w")
    fr = open("UNV/" + f + "/" + f + "_cap","r")
    for Ls in range(64):
        line = fr.readline()
        entries = line.split()
        for i in range(0,Ls):
            fw.write(entries[i] + ",")
        for i in range(Ls+1,64):
            fw.write(entries[i] + ",")
        fw.write("\n")
