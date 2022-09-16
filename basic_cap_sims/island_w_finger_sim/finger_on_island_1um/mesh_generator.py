
import sys
import os
import time
import platform
from datetime import datetime
import shutil
import salome
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import salome_notebook
import SMESH, SALOMEDS
from salome.smesh import smeshBuilder
import numpy as np

geompy = geomBuilder.New()

d_island = 0.05
delet = 0.03
cell = 3
d_air = 2
d_glass = 0.4
d_finger = .001;

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

## CREATE FACES

# bot layer

gen_bot = geompy.MakeFaceHW(cell,cell,1)
B0 = geompy.MakeTranslation(gen_bot,cell/2,cell/2,0)
Blist = [B0]

# fingerip

F0 = geompy.MakeTranslation(gen_bot,cell/2,cell/2,d_glass + d_finger)
Flist = [F0]

# top layer

gen_top = geompy.MakeFaceHW((cell - 2*delet - d_island)/2, cell, 1)
T0 = geompy.MakeTranslation(gen_top, (cell - 2*delet - d_island)/4, cell/2, d_glass)
T1 = geompy.MakeTranslation(T0, (cell + 2*delet + d_island)/2, 0, 0)
Tlist = [T0, T1];

# IX

gen_island = geompy.MakeFaceHW(d_island, cell, 1)
I = geompy.MakeTranslation(gen_island, cell/2, cell/2, d_glass)
Ilist = [I]

# create solids

box_air0 = geompy.MakeBoxDXDYDZ(cell,cell,d_air)
box_air_top = geompy.MakeTranslation(box_air0,0,0,d_glass+d_finger)
box_air_bot = geompy.MakeTranslation(box_air0,0,0,-d_air)
box_air1 = geompy.MakeBoxDXDYDZ(cell,cell,d_finger)
box_air_mid = geompy.MakeTranslation(box_air1,0,0,d_glass)

box_glass0 = geompy.MakeBoxDXDYDZ(cell,cell, d_glass)
box_glass = geompy.MakeTranslation(box_glass0,0,0,0)

# making the partitions

list_of_bodies = Tlist + Blist + Ilist + Tlist + [box_air_bot, box_air_top, box_air_mid, box_glass]

partition = geompy.MakePartition(list_of_bodies, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
[air_bot, glass, air_mid, air_top] = geompy.ExtractShapes(partition, geompy.ShapeType["SOLID"], True)
Faces = geompy.ExtractShapes(partition, geompy.ShapeType["FACE"], True)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

geompy.addToStudy(Flist[0], 'F0');
geompy.addToStudy(Tlist[0], 'T0');
geompy.addToStudy(Tlist[1], 'T1');
geompy.addToStudy(Ilist[0], 'I0');
geompy.addToStudy(Blist[0], 'B0');

geompy.addToStudy( box_air_bot, 'airb' )
geompy.addToStudy( box_air_top, 'airt' )
geompy.addToStudy( box_air_mid, 'airm' )
geompy.addToStudy( box_glass, 'glass' )
geompy.addToStudy( partition, 'partition' )

geompy.addToStudyInFather(partition,air_top,'airt');
geompy.addToStudyInFather(partition,air_bot,'airb');
geompy.addToStudyInFather(partition,air_mid,'airm');
geompy.addToStudyInFather(partition,glass,'glass');

# locating important Faces

Vz = geompy.MakeVectorDXDYDZ(0, 0, 1)
B_vertex = geompy.MakeVertex(0, 0, 0);
T_vertex = geompy.MakeVertex(0, 0, d_glass)
airt_vertex = geompy.MakeVertex(0, 0, d_glass+d_finger+d_air);
airb_vertex = geompy.MakeVertex(0, 0, -d_air)
airm_vertex = geompy.MakeVertex(0, 0, d_glass+d_finger)

face_group = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])

B_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, B_vertex, GEOM.ST_ON)
T_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, T_vertex, GEOM.ST_ON)
airb_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, airb_vertex, GEOM.ST_ON)
airt_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, airt_vertex, GEOM.ST_ON)
airm_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, airm_vertex, GEOM.ST_ON)

geompy.UnionList(face_group,T_plane + B_plane + airb_plane + airm_plane + airt_plane)

geompy.addToStudyInFather(partition,face_group,'face_group')
faces = geompy.ExtractShapes(face_group, geompy.ShapeType["FACE"], True)
print(len(faces))

face_names = ['left_C','airb','GND','island','finger','airt','right_C']
face_indexes = [0,2,3,4,5,6,8]
for i in range(len(face_indexes)):
    geompy.addToStudyInFather(face_group,faces[face_indexes[i]],face_names[i]);

# for i in range(len(faces)):
#     geompy.addToStudyInFather(face_group,faces[i],str(i))

#### create mesh groups

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(partition)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( .1 )
NETGEN_3D_Parameters_1.SetMinSize( .0001 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 4 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 208 )

NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom = face_group)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D_1.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( .1)
NETGEN_2D_Parameters_1.SetMinSize( .0001 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 4 )
NETGEN_2D_Parameters_1.SetChordalError( -1 )
NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_1.SetUseDelauney( 109 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_1.SetCheckChartBoundary( 192 )

isDone = Mesh_1.Compute()
Sub_mesh_1 = NETGEN_1D_2D_1.GetSubMesh()

airb = Mesh_1.GroupOnGeom(air_bot, 'air_bot', SMESH.VOLUME)
airt = Mesh_1.GroupOnGeom(air_top, 'air_top', SMESH.VOLUME)
airt = Mesh_1.GroupOnGeom(air_mid, 'air_mid', SMESH.VOLUME)
glass = Mesh_1.GroupOnGeom(glass, 'glass', SMESH.VOLUME)

all_faces = [];

for i in range(len(face_names)):
    all_faces.append(Mesh_1.GroupOnGeom(faces[face_indexes[i]], face_names[i], SMESH.FACE))

Mesh_1.ExportUNV('mesh.unv')

Groups = Mesh_1.GetGroups()

smesh.SetName(NETGEN_1D_2D_3D, 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_1D_2D_1, 'NETGEN 1D-2D-1')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
