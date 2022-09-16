#!/usr/bin/env python
# This file generates a series of FEM meshes which will then be used to compute capacitances for

import os
from os import listdir

hor = 1;
ver = 1;

indir = listdir("./")
dirs = []
for f in indir:
    if not os.path.isfile("./" + f):
        dirs.append(f)

for file in dirs:
    fw = open(file + "/" + file + ".sif","w")
    fw.write("""Header
  CHECK KEYWORDS Warn
  Mesh DB "." "."
  Include Path ""
  Results Directory ""
End

Simulation
  Max Output Level = 5
  Coordinate System = Cartesian
  Coordinate Mapping(3) = 1 2 3
  Simulation Type = Steady state
  Steady State Max Iterations = 1
  Output Intervals = 1
  Timestepping Method = BDF
  BDF Order = 1
  Coordinate Scaling = .001
  """
      )
    fw.write("Solver Input File = " + file + ".sif\n")
    fw.write("  Post File = " + file + ".vtu\n")
    fw.write("""End

Constants
  Gravity(4) = 0 -1 0 9.82
  Stefan Boltzmann = 5.67e-08
  Permittivity of Vacuum = 8.8542e-12
  Boltzmann Constant = 1.3807e-23
  Unit Charge = 1.602e-19
End

Body 1
  Target Bodies(1) = 1
  Name = "Body 1"
  Equation = 1
  Material = 2
End

Body 2
  Target Bodies(1) = 2
  Name = "Body 2"
  Equation = 1
  Material = 2
End

Body 3
  Target Bodies(1) = 3
  Name = "Body 3"
  Equation = 1
  Material = 4
End

Body 4
  Target Bodies(1) = 4
  Name = "Body 4"
  Equation = 1
  Material = 1
End

Body 5
  Target Bodies(1) = 5
  Name = "Body 5"
  Equation = 1
  Material = 3
End

Body 6
  Target Bodies(1) = 6
  Name = "Body 6"
  Equation = 1
  Material = 1
End

Solver 1
  Equation = Electrostatics
  """
          )
    fw.write("Capacitance Matrix Filename = " + file + "_cap" + "\n")
    fw.write(
    """  Variable = Potential
  Calculate Electric Flux = True
  Calculate Electric Field = True
  Calculate Electric Energy = True
  Procedure = "StatElecSolve" "StatElecSolver" """)
    fw.write("\n  Calculate Capacitance Matrix = True\n")
    fw.write("""  Exec Solver = Always
  Stabilize = True
  Bubbles = False
  Lumped Mass Matrix = False
  Optimize Bandwidth = True
  Steady State Convergence Tolerance = 1.0e-5
  Nonlinear System Convergence Tolerance = 1.0e-7
  Nonlinear System Max Iterations = 20
  Nonlinear System Newton After Iterations = 3
  Nonlinear System Newton After Tolerance = 1.0e-3
  Nonlinear System Relaxation Factor = 1
  Linear System Solver = Iterative
  Linear System Iterative Method = BiCGStab
  Linear System Max Iterations = 500
  Linear System Convergence Tolerance = 1.0e-10
  BiCGstabl polynomial degree = 2
  Linear System Preconditioning = ILU1
  Linear System ILUT Tolerance = 1.0e-3
  Linear System Abort Not Converged = False
  Linear System Residual Output = 10
  Linear System Precondition Recompute = 1
End

Equation 1
  Name = "Equation 1"
  Active Solvers(1) = 1
End

Material 1
  Name = "glass"
  Relative Permittivity = 7.3
End

Material 2
  Name = "Air"
  Relative Permittivity = 1.00059
End

Material 3
  Name = "PVB"
  Relative Permittivity = 2.5
End

Material 4
  Name = "OCA"
  Relative Permittivity = 4.0
End

""")

    fw.write("\n")
    fw.write("Boundary Condition 1\n")
    fw.write("  Target Boundaries(2) = 7 9\n")
    fw.write("  Name = \"air_ends\"\n")
    fw.write("  Electric Infinity BC = True\n")
    fw.write("End\n\n")

    fw.write("Boundary Condition 2\n")
    fw.write("  Target Boundaries(1) = 8\n")
    fw.write("  Name = \"GND\"\n")
    fw.write("  Capacitance Body = 0\n")
    fw.write("  Potential = 0\n")
    fw.write("End\n\n")

    bc = 3;
    bound = 10
    hn = 0;
    tn = 0;
    rn = 0;
    i_n = 0;
    iin = 0;
    cb = 1;

    for h in range(hor*3):

        fw.write("Boundary Condition " + str(bc) + "\n")
        fw.write("  Target Boundaries(1) = " + str(bound) + "\n")
        fw.write("  Name = \"Hx" + str(hn) + "\"\n")
        fw.write("  Capacitance Body = " + str(cb) + "\n")
        fw.write("End\n\n")

        bc += 1
        bound += 1
        hn += 1
        cb +=1

    bound_o = bound
    bound += (hor-1)*3 + 2

    iixc = [22,23,22]

    for h in range(hor*3):
        bounds = ""
        for i in range(iixc[h]*hor):
            bounds += " " + str(bound)
            bound += 1

        fw.write("Boundary Condition " + str(bc) + "\n")
        fw.write("  Target Boundaries(" + str(iixc[h]*ver) + ") = " + bounds + "\n")
        fw.write("  Name = \"IIx" + str(iin) + "\"\n")
        fw.write("  Capacitance Body = " + str(cb) + "\n")
        fw.write("End\n\n")
        bc += 1
        iin += 1
        cb +=1

    for h in range(hor*3):
        bounds = ""
        for i in range(6*ver):
            bounds += " " + str(bound)
            bound += 1

        fw.write("Boundary Condition " + str(bc) + "\n")
        fw.write("  Target Boundaries(" + str(ver*6) + ") =" + bounds + "\n")
        fw.write("  Name = \"Tx" + str(tn) + "\"\n")
        fw.write("  Capacitance Body = " + str(cb) + "\n")
        fw.write("End\n\n")

        bc += 1
        tn += 1
        cb +=1

    for v in range(ver*3):

        fw.write("Boundary Condition " + str(bc) + "\n")
        fw.write("  Target Boundaries(1) = " + str(bound) + "\n")
        fw.write("  Name = \"Rx" + str(rn) + "\"\n")
        fw.write("  Capacitance Body = " + str(cb) + "\n")
        fw.write("End\n\n")

        bc += 1
        bound += 1
        rn += 1
        cb +=1

    for h in range((hor-1)*3 + 2):

        fw.write("Boundary Condition " + str(bc) + "\n")
        fw.write("  Target Boundaries(1) = " + str(bound_o) + "\n")
        fw.write("  Name = \"Ix" + str(i_n) + "\"\n")
        fw.write("  Capacitance Body = " + str(cb) + "\n")
        fw.write("End\n\n")

        bc += 1
        bound_o += 1
        i_n += 1
        cb +=1

    fw.close()
