#!/usr/bin/env python
# This file generates a series of FEM meshes which will then be used to compute capacitances for

import os

unvdir = os.getcwd() + "/UNV/"
unvfiles = os.listdir(os.getcwd() + "/UNV/")

unvfiles = unvfiles[:-1:2]
for file in unvfiles:
    fr = open("sif_base.txt","r")
    fw = open(unvdir + file + "/" + file + ".sif","w")
    newline = fr.readline()
    while "*input_sif_file*" not in newline:
        print(newline)
        fw.write(newline)
        newline = fr.readline()
    fw.write("  Solver Input File = " + file + ".sif\n")
    fw.write("  Post File = " + file + ".vtu\n")
    newline = fr.readline()
    newline = fr.readline()
    while newline != "":
        fw.write(newline)
        newline = fr.readline()

    fr.close()

    # fr = open(unvdir + file + "/" + "entities.sif","r")
    # line = fr.readline()
    # print(line)
    # words = line.split()
    # while "boundaries" not in words:
    #     fr = open(unvdir + file + "/" + "entities.sif","r")
    #     line = fr.readline()
    #     words = line.split()
    #
    # while line != "":
    #     line = fr.readline()
    #     fw.wrie(line)
    #     line = fr.readline()
    #     fw.wrie(line)
    #     if "left_p" in line.split():
    #         fw.write("Potential = 1 \n")
    #         fw.write("Capacitance Body = 1 \n")
    #     elif "right_p" in line.split():
    #         fw.write("Capacitance Body = 2 \n")
    #         fw.write("Potential Constant = Logical True\n")
    #     else:
    #         fw.write("Potential Constant = Logical True\n")
    #
    #     line = fr.readline()
    #     fw.wrie(line)

    fw.write("Boundary Condition 1\n")
    fw.write("  Target Boundaries(11) = 1 2 4 5 6 8 9 10 12 13\n")
    fw.write("  Name = \"walls\"\n")
    fw.write("  Potential = 0\n")
    fw.write("End\n\n")

    fw.write("Boundary Condition 2\n")
    fw.write("  Target Boundaries(1) = 3\n")
    fw.write("  Name = \"left\"\n")
    fw.write("  Capacitance Body = 1\n")
    fw.write("  Potential = 1\n")
    fw.write("End\n\n")

    fw.write("Boundary Condition 2\n")
    fw.write("  Target Boundaries(1) = 11\n")
    fw.write("  Name = \"right\"\n")
    fw.write("  Capacitance Body = 2\n")
    fw.write("  Potential = -1\n")
    fw.write("End\n\n")

fw.close()
fr.close()
