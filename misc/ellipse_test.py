#!/usr/bin/env python

###
### This file is generated automatically by SALOME v9.7.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
import salome_notebook
notebook = salome_notebook.NoteBook()
sys.path.insert(0, r'/media/sf_Shared_Folder/misc')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New()

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Ellipse_1 = geompy.MakeEllipse(None, None, 200, 100)
Face_1 = geompy.MakeFaceWires([Ellipse_1], 1)
Rotation_1 = geompy.MakeRotation(Face_1, OZ, 67.5*math.pi/180.0)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Ellipse_1, 'Ellipse_1' )
geompy.addToStudy( Face_1, 'Face_1' )
geompy.addToStudy( Rotation_1, 'Rotation_1' )


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
