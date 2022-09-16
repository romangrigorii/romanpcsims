import numpy as np
import math
import os

unvdir = os.getcwd() + "/UNV/"
unvfiles = os.listdir(os.getcwd() + "/UNV/")
unvfiles = unvfiles[:-1:2]

for file in unvfiles:
    print(file)
    fr = open(unvdir + file + "/" + file + "_cap","r")
    print(fr.read())
