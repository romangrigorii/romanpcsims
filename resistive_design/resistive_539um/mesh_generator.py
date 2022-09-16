#!/usr/bin/env python

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

delet = 0.03
Hx_delet = 0.11
Rx_r = 1
Rx_d = 0.15

step = 5.39
cell_w_Hx = step - delet
cell_w_Tx = step - delet
cell_h = step/2

d_air = 10
d_glass = .4
d_PVB = 0.33
d_spine = 0.7
d_OCA = 0.25

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)

ver = 1
hor = 1

## CREATE FACES

# HX

gen_cell_Hx = geompy.MakeFaceHW(cell_w_Hx, step*ver, 1)
Hx0t = geompy.MakeTranslation(gen_cell_Hx, step/2, step/2, d_glass)
Hx0 = Hx0t
if Hxtype == 1:
    deletion_line = geompy.MakeFaceHW(cell_w_Hx - Hxthick, Hxdel, 1)
    sqcut_list = [];
    for y in range(0,int(step*ver/(Hxdel+Hxthick))):
        if np.mod(y,2) == 0:
            sqcut_list.append(geompy.MakeTranslation(deletion_line,step/2 + Hxthick/2,(y+1)*(Hxthick+Hxdel),d_glass))
        else:
            sqcut_list.append(geompy.MakeTranslation(deletion_line,step/2 - Hxthick/2,(y+1)*(Hxthick+Hxdel),d_glass))

    Hx0 = geompy.MakeCutList(Hx0t, sqcut_list, True)
Hx1 = geompy.MakeTranslation(Hx0, step, 0,0)
Hxlist = [Hx0,Hx1]

for h in range(hor-1):
    Hxlist.append(geompy.MakeTranslation(Hxlist[h*2], step, 0, 0))
    Hxlist.append(geompy.MakeTranslation(Hxlist[h*2 + 1], step, 0, 0))


# RX

gen_rx_c = geompy.MakeDiskR(Rx_r, 1)
gen_rx_b = geompy.MakeFaceHW(step, Rx_d, 1)
gen_rx_cb = geompy.MakeFuseList([gen_rx_c, gen_rx_b], True, True)

Rx0_0 = geompy.MakeTranslation(gen_rx_cb,step/2, step/4, 0)
Rx0_1 = geompy.MakeTranslation(gen_rx_b,step/2, step*3/4, 0)
Rx1_0 = geompy.MakeTranslation(gen_rx_b,step*3/2, step/4, 0)
Rx1_1 = geompy.MakeTranslation(gen_rx_cb,step*3/2, step*3/4, 0)

Rx0 = geompy.MakeFuseList([Rx0_0, Rx1_0], True, True)
Rx1 = geompy.MakeFuseList([Rx0_1, Rx1_1], True, True)

Rxlist0 = [Rx0]
Rxlist1 = [Rx1]

for h in range(hor-1):
    Rxlist0.append(geompy.MakeTranslation(Rxlist0[-1], step*2, 0, 0))
    Rxlist1.append(geompy.MakeTranslation(Rxlist1[-1], step*2, 0, 0))

Rx0 = geompy.MakeFuseList(Rxlist0, True, True)
Rx1 = geompy.MakeFuseList(Rxlist1, True, True)

Rxlist = [Rx0, Rx1]

for v in range(ver-1):
    Rxlist.append(geompy.MakeTranslation(Rxlist[v*2], 0, step, 0))
    Rxlist.append(geompy.MakeTranslation(Rxlist[v*2 + 1], 0, step, 0))

# TX

gen_rx_deletion_c = geompy.MakeDiskR(Rx_r+delet,1)
gen_rx_deletion_b = geompy.MakeFaceHW(step, Rx_d + delet*2, 1)
gen_rx_deletion_cb = geompy.MakeFuseList([gen_rx_deletion_c,gen_rx_deletion_b])

gen_cell_Tx = geompy.MakeFaceHW(cell_w_Tx, cell_h, 1)

gen_tx_0 = geompy.MakeCutList(gen_cell_Tx, [gen_rx_deletion_cb], True)
gen_tx_1 = geompy.MakeCutList(gen_cell_Tx, [gen_rx_deletion_b], True)

Tx0_0 = geompy.MakeTranslation(gen_tx_0,step/2, step/4, 0)
Tx0_1 = geompy.MakeTranslation(gen_tx_1,step/2, step*3/4, 0)
Tx1_0 = geompy.MakeTranslation(gen_tx_1,step*3/2, step/4, 0)
Tx1_1 = geompy.MakeTranslation(gen_tx_0,step*3/2, step*3/4, 0)

Tx0t = geompy.MakeFuseList([Tx0_0, Tx0_1], True, True)
Tx1t = geompy.MakeFuseList([Tx1_0, Tx1_1], True, True)

Txlist0 = [Tx0t]
Txlist1 = [Tx1t]

for v in range(ver-1):
    Txlist0.append(geompy.MakeTranslation(Txlist0[v], 0, step, 0))
    Txlist1.append(geompy.MakeTranslation(Txlist1[v], 0, step, 0))

Tx0 = geompy.MakeFuseList(Txlist0, True, True)
Tx1 = geompy.MakeFuseList(Txlist1, True, True)

Txlist = [Tx0, Tx1]

for h in range(hor-1):
    Txlist.append(geompy.MakeTranslation(Txlist[h*2], step*2, 0, 0))
    Txlist.append(geompy.MakeTranslation(Txlist[h*2 + 1], step*2, 0, 0))

# create solids

box_air0 = geompy.MakeBoxDXDYDZ(step*2*hor,step*ver,d_air+d_glass+d_PVB+d_spine+d_OCA)
box_air1 = geompy.MakeTranslation(box_air0,0,0, -(d_air + d_glass + d_PVB + d_spine + d_OCA)/2 - ((d_PVB + d_spine + d_OCA - d_glass)/2))

box_glass0 = geompy.MakeBoxDXDYDZ(step*2*hor,step*ver,d_glass)
box_glass = geompy.MakeTranslation(box_glass0,0,0,0)

box_PVB0 = geompy.MakeBoxDXDYDZ(step*2*hor,step*ver,d_PVB)
box_PVB = geompy.MakeTranslation(box_PVB0,0,0,-d_PVB)

box_spine0 = geompy.MakeBoxDXDYDZ(step*2*hor,step*ver,d_spine)
box_spine = geompy.MakeTranslation(box_spine0,0,0,-d_PVB - d_spine)

box_OCA0 = geompy.MakeBoxDXDYDZ(step*2*hor,step*ver,d_OCA)
box_OCA = geompy.MakeTranslation(box_OCA0,0,0,-d_PVB - d_spine - d_OCA)

box_air = geompy.MakeCutList(box_air1,[box_glass,box_PVB,box_spine,box_OCA],True)

# making the partitions

list_of_bodies = Rxlist + Txlist + Hxlist + [box_OCA,box_PVB,box_air,box_glass,box_spine]

partition = geompy.MakePartition(list_of_bodies, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
[airb,OCA,spine,PVB,glass,airt] = geompy.ExtractShapes(partition, geompy.ShapeType["SOLID"], True)
Faces = geompy.ExtractShapes(partition, geompy.ShapeType["FACE"], True)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

for h in range(hor):
    geompy.addToStudy( Hxlist[h*2], 'Hx' + str(h*2))
    geompy.addToStudy( Hxlist[h*2 + 1], 'Hx' + str(h*2 + 1))

for v in range(ver):
    geompy.addToStudy( Rxlist[v*2], 'Rx' + str(v*2))
    geompy.addToStudy( Rxlist[v*2 + 1], 'Rx' + str(v*2 + 1))

for h in range(hor):
    geompy.addToStudy( Txlist[h*2], 'Tx' + str(h*2))
    geompy.addToStudy( Txlist[h*2 + 1], 'Tx' + str(h*2 + 1))

geompy.addToStudy( box_air, 'air' )
geompy.addToStudy( box_glass, 'glass' )
geompy.addToStudy( box_PVB, 'PVB' )
geompy.addToStudy( box_spine, 'spine' )
geompy.addToStudy( box_OCA, 'OCA' )
geompy.addToStudy( partition, 'partition' )

geompy.addToStudyInFather(partition,airt,'airt');
geompy.addToStudyInFather(partition,glass,'glass');
geompy.addToStudyInFather(partition,PVB,'PVB');
geompy.addToStudyInFather(partition,spine,'spine');
geompy.addToStudyInFather(partition,OCA,'OCA');
geompy.addToStudyInFather(partition,airb,'airb');

# locating important Faces

Vz = geompy.MakeVectorDXDYDZ(0, 0, 1)
GND_vertex = geompy.MakeVertex(0, 0, -(d_PVB+d_spine+d_OCA))
Tx_vertex = geompy.MakeVertex(0, 0, 0)
Hx_vertex = geompy.MakeVertex(0, 0, d_glass)
airt_vertex = geompy.MakeVertex(0, 0, d_air/2 + d_glass);
airb_vertex = geompy.MakeVertex(0, 0, -(d_air/2 + d_PVB + d_spine + d_OCA))

face_group0 = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])
face_group1 = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])
face_group2 = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])

GND_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, GND_vertex, GEOM.ST_ON)
airt_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, airt_vertex, GEOM.ST_ON)
airb_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, airb_vertex, GEOM.ST_ON)
Hx_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, Hx_vertex, GEOM.ST_ON)
Tx_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                    Vz, Tx_vertex, GEOM.ST_ON)

geompy.UnionList(face_group0,[airb_plane[0],GND_plane[0],airt_plane[0]])
geompy.UnionList(face_group1,Hx_plane)
geompy.UnionList(face_group2,Tx_plane)

geompy.addToStudyInFather(partition,face_group0,'face_group0')
faces0 = geompy.ExtractShapes(face_group0, geompy.ShapeType["FACE"], True)
face0_names = ['airb_face','GND_face','airt_face']
for i in range(len(face0_names)):
    geompy.addToStudyInFather(face_group0,faces0[i],face0_names[i])

# identifying the faces and grouping, Hx layer

geompy.addToStudyInFather(partition,face_group1,'face_group1')
faces1 = geompy.ExtractShapes(face_group1, geompy.ShapeType["FACE"], True)

face1_names_Hx = []
face1_indexes_Hx = []

for i in range(hor*2):
    face1_names_Hx.append('Hx' + str(i))
for i in range(hor*2):
    face1_indexes_Hx.append(i*2+1)

face1_names = face1_names_Hx
face1_indexes = face1_indexes_Hx
faces1_new = [];

for i in face1_indexes:
    faces1_new.append(faces1[i])

for i in range(len(faces1_new)):
    geompy.addToStudyInFather(face_group1,faces1_new[i],face1_names[i])

# identifying the faces and grouping, Tx layer

geompy.addToStudyInFather(partition,face_group2,'face_group2')
faces2 = geompy.ExtractShapes(face_group2, geompy.ShapeType["FACE"], True)

face2_names_Tx = [];
face2_indexes_Tx = [];
face2_names_Rx = [];
face2_indexes_Rx = [];

for h in range(hor*2):
    for v in range(ver*2*2):
        face2_names_Tx.append('Tx' + str(h) + '_' + str(v))
for v in range(ver):
    face2_names_Rx.append('Rx' + str(v*2))
for v in range(ver):
    face2_names_Rx.append('Rx' + str(v*2 + 1))

if ver == 1:
    L1 = hor*ver*4
    L2 = L1 + 4
    L3 = L2 + 1
    L4 = L3 + L1
    for i in range(L1):
        face2_indexes_Tx.append(i)
    for i in range(L3,L4):
        face2_indexes_Tx.append(i)
    face2_indexes_Rx = [L1,L2]
else:
    L1 = ver*4*(hor-1) + 6
    L5 = ver*4*hor*2 + (ver*4) + 1
    L4 = L5 - L1
    L2 = int((L5 - 2*ver - 3)/2)
    L3 = int((L5 + 2*ver + 1)/2)
    # print([L1,L2,L3,L4,L5])

    for i in range(L1):
        face2_indexes_Tx.append(i)

    for i in range(L1,L2):
        if np.mod(i - L1,5) == 0:
            face2_indexes_Rx.append(i)
        else:
            face2_indexes_Tx.append(i)

    face2_indexes_Rx.append(L2)
    face2_indexes_Rx.append(L3)

    faceindRx = []
    faceindTx = []

    for i in range(1,L1):
        faceindTx.append(L5 - i)

    for i in range(L1,L2+1):
        if np.mod(i - L1 - 1,5) == 0:
            faceindRx.append(L5 - i)
        else:
            faceindTx.append(L5 - i)

    faceindTx.reverse()
    faceindRx.reverse()

    face2_indexes_Rx += faceindRx
    face2_indexes_Tx += faceindTx

face2_names = face2_names_Tx + face2_names_Rx
face2_indexes = face2_indexes_Tx + face2_indexes_Rx
print(face2_indexes)
print(len(faces2))
faces2_new = [];
for i in face2_indexes:
    faces2_new.append(faces2[i])

for i in range(len(face2_names)):
    geompy.addToStudyInFather(face_group2,faces2_new[i],face2_names[i])

for i in range(len(faces2)):
    geompy.addToStudyInFather(face_group2,faces2[i],'A_' + str(i))

# create mesh groups

smesh = smeshBuilder.New()
Mesh_1 = smesh.Mesh(partition)
NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
NETGEN_3D_Parameters_1.SetMaxSize( 1 )
NETGEN_3D_Parameters_1.SetMinSize( .01 )
NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
NETGEN_3D_Parameters_1.SetOptimize( 1 )
NETGEN_3D_Parameters_1.SetFineness( 4 )
NETGEN_3D_Parameters_1.SetChordalError( -1 )
NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
NETGEN_3D_Parameters_1.SetCheckChartBoundary( 208 )

NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=face_group0)
NETGEN_2D_Parameters_1 = NETGEN_1D_2D_1.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( .1)
NETGEN_2D_Parameters_1.SetMinSize( .001 )
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

NETGEN_1D_2D_2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=face_group1)
NETGEN_2D_Parameters_2 = NETGEN_1D_2D_2.Parameters()
NETGEN_2D_Parameters_2.SetMaxSize( .1)
NETGEN_2D_Parameters_2.SetMinSize( .001 )
NETGEN_2D_Parameters_2.SetSecondOrder( 0 )
NETGEN_2D_Parameters_2.SetOptimize( 1 )
NETGEN_2D_Parameters_2.SetFineness( 4 )
NETGEN_2D_Parameters_2.SetChordalError( -1 )
NETGEN_2D_Parameters_2.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_2.SetFuseEdges( 1 )
NETGEN_2D_Parameters_2.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_2.SetUseDelauney( 109 )
NETGEN_2D_Parameters_2.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_2.SetCheckChartBoundary( 192 )

NETGEN_1D_2D_3 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=face_group2)
NETGEN_2D_Parameters_3 = NETGEN_1D_2D_3.Parameters()
NETGEN_2D_Parameters_3.SetMaxSize( .1 )
NETGEN_2D_Parameters_3.SetMinSize( .001 )
NETGEN_2D_Parameters_3.SetSecondOrder( 0 )
NETGEN_2D_Parameters_3.SetOptimize( 1 )
NETGEN_2D_Parameters_3.SetFineness( 4 )
NETGEN_2D_Parameters_3.SetChordalError( -1 )
NETGEN_2D_Parameters_3.SetChordalErrorEnabled( 0 )
NETGEN_2D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_3.SetFuseEdges( 1 )
NETGEN_2D_Parameters_3.SetWorstElemMeasure( 0 )
NETGEN_2D_Parameters_3.SetUseDelauney( 109 )
NETGEN_2D_Parameters_3.SetQuadAllowed( 0 )
NETGEN_2D_Parameters_3.SetCheckChartBoundary( 192 )

isDone = Mesh_1.Compute()
Sub_mesh_1 = NETGEN_1D_2D_1.GetSubMesh()
Sub_mesh_2 = NETGEN_1D_2D_2.GetSubMesh()
Sub_mesh_3 = NETGEN_1D_2D_3.GetSubMesh()

#
airb1 = Mesh_1.GroupOnGeom(airb, 'airb', SMESH.VOLUME)
airt1 = Mesh_1.GroupOnGeom(airt, 'airt', SMESH.VOLUME)
OCA1 = Mesh_1.GroupOnGeom(OCA, 'OCA', SMESH.VOLUME)
spine1 = Mesh_1.GroupOnGeom(spine, 'spine', SMESH.VOLUME)
PVB1 = Mesh_1.GroupOnGeom(PVB, 'PVB', SMESH.VOLUME)
glass1 = Mesh_1.GroupOnGeom(glass, 'glass', SMESH.VOLUME)

all_faces = [];

for i in range(len(face0_names)):
    all_faces.append(Mesh_1.GroupOnGeom(faces0[i], face0_names[i], SMESH.FACE))
for i in range(len(face1_names)):
    all_faces.append(Mesh_1.GroupOnGeom(faces1_new[i], face1_names[i], SMESH.FACE))
for i in range(len(face2_names)):
    all_faces.append(Mesh_1.GroupOnGeom(faces2_new[i], face2_names[i], SMESH.FACE))

Mesh_1.ExportUNV('mesh.unv')

Groups = Mesh_1.GetGroups()

smesh.SetName(NETGEN_1D_2D_3D, 'NETGEN 1D-2D-3D')
smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_1D_2D_1, 'NETGEN 1D-2D-1')
smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 3D Parameters_1')
smesh.SetName(NETGEN_1D_2D_2, 'NETGEN 1D-2D-2')
smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 3D Parameters_2')
smesh.SetName(NETGEN_1D_2D_3, 'NETGEN 1D-2D-3')
smesh.SetName(NETGEN_2D_Parameters_3, 'NETGEN 3D Parameters_3')

for i in range(6):
    smesh.SetName(Groups[i],'Solid_' + str(i))
for i in range(6,len(Groups)):
    smesh.SetName(Groups[i],'Face_' + str(i-6))

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
