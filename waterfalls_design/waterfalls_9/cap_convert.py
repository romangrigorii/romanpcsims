# this file converts the mesh_cap to maxwell_cap
import csv
import numpy as np
import os
from os import listdir
import math

indir = listdir("./")
dirs = []
for f in indir:
    if not os.path.isfile("./" + f):
        dirs.append(f)

for fname in dirs:
    print(fname)
    file = open(fname + '/' + fname + '_cap','r')
    line = file.readline();
    nums = line.split();
    L = len(nums) + 1
    mat = np.zeros((L,L))
    v = 0
    f = open(fname + '/mesh_maxwell', 'w')
    writer = csv.writer(f)

    while line != "":
        nums = line.split();
        print(nums)
        for h in range(L-1):
            if h == v:
                mat[h,L-1] = -float(nums[h])/1e-12
                mat[L-1,h] = -float(nums[h])/1e-12
            else:
                mat[v,h] = -float(nums[h])/1e-12
        line = file.readline();
        v+=1

    for v in range(L):
        for h in range(L):
            if h!=v:
                mat[v,v] -= mat[v,h]
        writer.writerow(mat[v,:])
    f.close()
