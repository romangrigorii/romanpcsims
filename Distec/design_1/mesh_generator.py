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

xlen = 3
ylen = 3

w_island = 0.4 # width of the island
L_island = 4.495 # length of the islands
delet = 0.03 # deletion line
delet2 = .05 # deletion line between Tx
Rx_L = 2.0075 # length of the rx head
Rx_d = .58 # width of the rx head
Rx_trunk = .27 # width of the Rx Trunk

step_w = 4.645
step_h = 4.645

cell_w_Hx = 1.7425
cell_w_Tx = 1.6625
cell_h_Hx = step_h
cell_h_Tx = step_h

d_air = 5
d_sensor = .4
d_OCA = .25
d_spine = 0.7
d_PVB = 0.38

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)

## CREATE FACES

# HX : these are Hx faces

gen_step_wx = geompy.MakeFaceHW(cell_w_Hx, step_h*ylen, 1)
gen_step_wx = geompy.MakeFuseList([geompy.MakeTranslation(gen_step_wx,cell_w_Hx/2,step_h*ylen/2,d_sensor), geompy.MakeTranslation(gen_step_wx,cell_w_Hx/2 + step_w/2,step_h*ylen/2,d_sensor)])
Hxlist = [gen_step_wx]
for i in range(xlen-1):
    Hxlist.append(geompy.MakeTranslation(Hxlist[-1], step_w, 0, 0))

# IX these are islands around the IIx

base1 = geompy.MakeFaceHW(w_island + 4*delet,L_island + 4*delet,1)
base2 = geompy.MakeFaceHW(w_island + 2*delet,L_island + 2*delet,1)
base = geompy.MakeCutList(base1,[base2])
base = geompy.MakeFuseList([geompy.MakeTranslation(base,(w_island + 6*delet)/2 + cell_w_Hx,cell_h_Hx/2,d_sensor),geompy.MakeTranslation(base,(w_island + 6*delet)/2 + cell_w_Hx + step_w/2 ,cell_h_Hx/2,d_sensor)],True,True)
Ixlist = []
for x in range(xlen):
    for y in range(ylen):
        Ixlist.append(geompy.MakeTranslation(base,x*step_w,y*step_h,0))

# TXIX
# none
# IIX

base = geompy.MakeFaceHW(w_island,L_island,1)
IIx = geompy.MakeFuseList([geompy.MakeTranslation(base,(w_island + 6*delet)/2 + cell_w_Hx,cell_h_Hx/2,d_sensor),geompy.MakeTranslation(base,(w_island + 6*delet)/2 + cell_w_Hx + step_w/2,cell_h_Hx/2,d_sensor)],True,True)
IIxlist = []
for x in range(xlen):
    for y in range(ylen):
        IIxlist.append(geompy.MakeTranslation(IIx,x*step_w,y*step_h,0))
print(len(IIxlist))

# ## creating Rx

gen_rx_t1 = geompy.MakeTranslation(geompy.MakeFaceHW(step_w, Rx_trunk, 1),step_w/2,(Rx_trunk + delet)/2,0)
gen_rx_t2 = geompy.MakeTranslation(geompy.MakeFaceHW(step_w, Rx_trunk, 1),step_w/2,-(Rx_trunk + delet)/2,0)
gen_rx_f = geompy.MakeFaceHW(Rx_d,Rx_L,1)
gen_rx_f1 = geompy.MakeFuseList([geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx,(Rx_trunk + delet/2 + Rx_L/2),0),geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx + step_w/2,(Rx_trunk + delet/2 + Rx_L/2),0)],True,True);
gen_rx_f2 = geompy.MakeFuseList([geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx,-(Rx_trunk + delet/2 + Rx_L/2),0),geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx + step_w/2,-(Rx_trunk + delet/2 + Rx_L/2),0)], True, True);
gen_rx_ft_1 = geompy.MakeFuseList([gen_rx_t1, gen_rx_f1], True, True)
gen_rx_ft_2 = geompy.MakeFuseList([gen_rx_t2, gen_rx_f2], True, True)
Rx1 = geompy.MakeFuseList([gen_rx_ft_1,geompy.MakeTranslation(gen_rx_ft_1,step_w,0,0),geompy.MakeTranslation(gen_rx_t1,2*step_w,0,0),geompy.MakeTranslation(gen_rx_t1,3*step_w,0,0)],True,True)
Rx2 = geompy.MakeFuseList([gen_rx_t2,geompy.MakeTranslation(gen_rx_t2,step_w,0,0),geompy.MakeTranslation(gen_rx_ft_2,2*step_w,0,0),geompy.MakeTranslation(gen_rx_ft_2,3*step_w,0,0)],True,True)
Rx0 = geompy.MakeFuseList([Rx1,Rx2],True,True)
Rxlist = [geompy.MakeTranslation(Rx0,0,step_h/2,0)]
for i in range(ylen-1):
    Rxlist.append(geompy.MakeTranslation(Rxlist[-1], 0, step_h, 0))

#
# # TX
#

gen_rx_t1 = geompy.MakeTranslation(geompy.MakeFaceHW(step_w, Rx_trunk+delet*2, 1),step_w/2,(Rx_trunk + delet)/2,0)
gen_rx_t2 = geompy.MakeTranslation(geompy.MakeFaceHW(step_w, Rx_trunk+delet*2, 1),step_w/2,-(Rx_trunk + delet)/2,0)
gen_rx_f = geompy.MakeFaceHW(Rx_d+delet*2,Rx_L+delet*2,1)
gen_rx_f1 = geompy.MakeFuseList([geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx,(Rx_trunk + delet + Rx_L/2),0),geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx + step_w/2,(Rx_trunk + delet + Rx_L/2),0)],True,True);
gen_rx_f2 = geompy.MakeFuseList([geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx,-(Rx_trunk + delet + Rx_L/2),0),geompy.MakeTranslation(gen_rx_f,(w_island + 6*delet)/2 + cell_w_Hx + step_w/2,-(Rx_trunk + delet + Rx_L/2),0)],True,True);
gen_rx_ft_1 = geompy.MakeFuseList([gen_rx_t1, gen_rx_f1], True, True)
gen_rx_ft_2 = geompy.MakeFuseList([gen_rx_t2, gen_rx_f2], True, True)
Rx1 = geompy.MakeFuseList([gen_rx_ft_1,geompy.MakeTranslation(gen_rx_ft_1,step_w,0,0),geompy.MakeTranslation(gen_rx_t1,2*step_w,0,0),geompy.MakeTranslation(gen_rx_t1,3*step_w,0,0)],True,True)
Rx2 = geompy.MakeFuseList([gen_rx_t2,geompy.MakeTranslation(gen_rx_t2,step_w,0,0),geompy.MakeTranslation(gen_rx_ft_2,2*step_w,0,0),geompy.MakeTranslation(gen_rx_ft_2,3*step_w,0,0)],True,True)
Rx0 = geompy.MakeFuseList([Rx1,Rx2],True,True)
RxlistTx = [geompy.MakeTranslation(Rx0,0,step_h/2,0)]
for i in range(ylen-1):
    RxlistTx.append(geompy.MakeTranslation(RxlistTx[-1], 0, step_h, 0))


Tx = geompy.MakeTranslation(geompy.MakeFaceHW(step_w - delet2, step_h*ylen, 1),step_w/2,step_h*ylen/2,0)
Txlist = [Tx]
for i in range(xlen-1):
    Txlist.append(geompy.MakeTranslation(Txlist[-1], step_w, 0, 0))
Tx = geompy.MakeFuseList(Txlist, True, True)
Txlist = [geompy.MakeCutList(Tx,RxlistTx)]


# ### create solids
#

box_air_surr = geompy.MakeBoxDXDYDZ(step_w*xlen,step_h*ylen,d_air)
box_air1 = geompy.MakeTranslation(box_air_surr, 0, 0, d_sensor)
box_air2 = geompy.MakeTranslation(box_air_surr, 0, 0, -(d_PVB+d_spine+d_OCA+d_air))

box_sensor = geompy.MakeBoxDXDYDZ(step_w*xlen,step_h*ylen,d_sensor)

box_PVB0 = geompy.MakeBoxDXDYDZ(step_w*xlen,step_h*ylen,d_PVB)
box_PVB = geompy.MakeTranslation(box_PVB0,0,0,-d_PVB)

box_spine0 = geompy.MakeBoxDXDYDZ(step_w*xlen,step_h*ylen,d_spine)
box_spine = geompy.MakeTranslation(box_spine0,0,0,-d_PVB-d_spine)

box_OCA0 = geompy.MakeBoxDXDYDZ(step_w*xlen,step_h*ylen,d_OCA)
box_OCA = geompy.MakeTranslation(box_OCA0,0,0,-d_PVB - d_spine - d_OCA)

# # making the partitions
list_of_bodies = Ixlist + Hxlist + IIxlist + Rxlist + Txlist + [box_OCA,box_PVB,box_air1,box_air2,box_sensor,box_spine]
partition = geompy.MakePartition(list_of_bodies, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
[airb,OCA,spine,PVB,sensor,airt] = geompy.ExtractShapes(partition, geompy.ShapeType["SOLID"], True)
Faces = geompy.ExtractShapes(partition, geompy.ShapeType["FACE"], True)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )


for i in range(len(Hxlist)):
    geompy.addToStudy( Hxlist[i], 'Hx' + str(i))
for i in range(len(IIxlist)):
    geompy.addToStudy( IIxlist[i], 'IIx' + str(i))
for i in range(len(Ixlist)):
    geompy.addToStudy( Ixlist[i], 'Ix' + str(i))
for i in range(len(Rxlist)):
    geompy.addToStudy( Rxlist[i], 'Rx' + str(i))
for i in range(len(Txlist)):
    geompy.addToStudy( Txlist[i], 'Tx' + str(i))


geompy.addToStudy( box_air1, 'airt' )
geompy.addToStudy( box_air2, 'airb' )
geompy.addToStudy( box_sensor, 'sensor' )
geompy.addToStudy( box_PVB, 'PVB' )
geompy.addToStudy( box_spine, 'spine' )
geompy.addToStudy( box_OCA, 'OCA' )
geompy.addToStudy( partition, 'partition' )

geompy.addToStudyInFather(partition,airt,'airt');
geompy.addToStudyInFather(partition,sensor,'sensor');
geompy.addToStudyInFather(partition,PVB,'PVB');
geompy.addToStudyInFather(partition,spine,'spine');
geompy.addToStudyInFather(partition,OCA,'OCA');
geompy.addToStudyInFather(partition,airb,'airb');
#
# # locating important Faces

Vz = geompy.MakeVectorDXDYDZ(0, 0, 1)
GND_vertex = geompy.MakeVertex(0, 0, -(d_PVB+d_spine+d_OCA))
Tx_vertex = geompy.MakeVertex(0, 0, 0)
Hx_vertex = geompy.MakeVertex(0, 0, d_sensor)
airt_vertex = geompy.MakeVertex(0, 0, d_air + d_sensor);
airb_vertex = geompy.MakeVertex(0, 0, -(d_air + d_PVB + d_spine + d_OCA))

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
#
# #### identifying the faces and grouping, Hx layer
#
geompy.addToStudyInFather(partition,face_group1,'face_group1')
faces1 = geompy.ExtractShapes(face_group1, geompy.ShapeType["FACE"], True)
print(len(faces1))

# face1_names_Hx = [];
# face1_indexes_Hx = [];
# face1_names_Ix = [];
# face1_indexes_Ix = [];
# face1_names_IIx = [];
# face1_indexes_IIx = [];
# for i in range(len(faces1)):
#     face1_names_Hx.append(str(i))
#     face1_indexes_Hx.append(i)

# face1_names_Hx = [];
# face1_indexes_Hx = [0,14,28,42,56,70,84,98];
# face1_names_Ix = [];
# face1_indexes_Ix = [2,5,9,12];
# face1_names_IIx = [];
# face1_indexes_IIx = [3,6,10,13];
# for i in range(len(face1_indexes_Hx)):
#     face1_names_Hx.append('Hx' + str(i))
# for i in range(len(face1_indexes_Ix)):
#     face1_names_Ix.append('Ix' + str(i))
# for i in range(len(face1_indexes_IIx)):
#     face1_names_IIx.append('IIx' + str(i))

face1_names_Hx = [];
face1_indexes_Hx = [];
face1_names_Ix = [];
face1_indexes_Ix = [];
face1_names_IIx = [];
face1_indexes_IIx = [];
hxi = 0;
iixi = 0;
ixi = 0;

for i in range(len(faces1)):
    [xmi,xma,ymi,yma,zmi,zma] = geompy.BoundingBox(faces1[i])
    if (np.abs(xmi-xma)< w_island*1.02) and (np.abs(xmi-xma) > w_island*.98):
            face1_indexes_Ix.append(i)
            face1_names_Ix.append('Ix' + str(ixi))
            ixi+=1
    if (np.abs(xmi-xma)< (w_island + 4*delet)*1.02) and (np.abs(xmi-xma) > (w_island + 4*delet)*.95):
        if (np.abs(ymi - yma) < (L_island + 4*delet)*1.02):
            face1_indexes_IIx.append(i)
            face1_names_IIx.append('IIx' + str(iixi))
            iixi+=1
    if (np.abs(xmi - xma) <  cell_w_Hx*1.02) and (np.abs(xmi - xma) > cell_w_Hx*.98):
        face1_indexes_Hx.append(i)
        face1_names_Hx.append('Hx' + str(hxi))
        hxi+=1


face1_names = face1_names_Hx + face1_names_Ix + face1_names_IIx
face1_indexes = face1_indexes_Hx + face1_indexes_Ix + face1_indexes_IIx
faces1_new = [];
for i in face1_indexes:
    faces1_new.append(faces1[i])
for i in range(len(faces1_new)):
    geompy.addToStudyInFather(face_group1,faces1_new[i],face1_names[i])
#
# #### identifying the faces and grouping, Tx layer
#
geompy.addToStudyInFather(partition,face_group2,'face_group2')
faces2 = geompy.ExtractShapes(face_group2, geompy.ShapeType["FACE"], True)
#
# face2_names_TxIx = [];
# face2_indexes_TxIx = [];
# face2_names_Tx = [];
# face2_indexes_Tx = [];
# face2_names_Rx = [];
# face2_indexes_Rx = [];
# for i in range(len(faces2)):
#     face2_names_Tx.append(str(i))
#     face2_indexes_Tx.append(i)

# face2_names_TxIx = [];
# face2_indexes_TxIx = [20,21,25,28,30]
# face2_names_Tx = [];
# face2_indexes_Tx = [1,4,8,9,12,13,14,15,16,17,23,26, 24,27,34,35,36,37,38,39,42,43,44,45]
# face2_names_Rx = [];
# face2_indexes_Rx = [11,33,18,40];
# for i in range(len(face2_indexes_Tx)):
#     face2_names_Tx.append("Tx" + str(i))
# for i in range(len(face2_indexes_Rx)):
#     face2_names_Rx.append("Rx" + str(i))
# for i in range(len(face2_indexes_TxIx)):
#     face2_names_TxIx.append("TxIx" + str(i))
#

face2_names_TxIx = [];
face2_indexes_TxIx = [];
face2_names_Tx = [];
face2_indexes_Tx = [];
face2_names_Rx = [];
face2_indexes_Rx = [];
txi = 0;
rxi = 0;
tixi = 0;
for i in range(len(faces2)):
    [xmi,xma,ymi,yma,zmi,zma] = geompy.BoundingBox(faces2[i])
    # if (np.abs(xmi-xma)< w_island*1.05) and (np.abs(xmi-xma) > w_island*.95):
    #     face2_indexes_TxIx.append(i)
    #     face2_names_TxIx.append('TxIx' + str(tixi))
    #     tixi+=1
    if (np.abs(xmi - xma) < step_w*1.05) and (np.abs(xmi - xma) > step_w/5):
        if (np.abs(xmi - xma) < step_h*1.05) and (np.abs(xmi - xma) > step_h/5):
            face2_indexes_Tx.append(i)
            face2_names_Tx.append('Tx' + str(txi))
            txi+=1
    if (np.abs(xmi - xma) <  step_w*xlen*1.05) and (np.abs(xmi - xma) > step_w*xlen*.95):
        if (np.abs(ymi - yma) < step_h/2*1.05) and (np.abs(ymi - yma) > step_h/2*.9):
            face2_indexes_Rx.append(i)
            face2_names_Rx.append('Rx' + str(rxi))
            rxi+=1

face2_names = face2_names_Tx + face2_names_Rx + face2_names_TxIx
face2_indexes = face2_indexes_Tx + face2_indexes_Rx + face2_indexes_TxIx
faces2_new = [];
for i in face2_indexes:
    faces2_new.append(faces2[i])

for i in range(len(face2_names)):
    geompy.addToStudyInFather(face_group2,faces2_new[i],face2_names[i])

## create mesh groups
#
# smesh = smeshBuilder.New()
# Mesh_1 = smesh.Mesh(partition)
# NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
# NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
# NETGEN_3D_Parameters_1.SetMaxSize( 1)
# NETGEN_3D_Parameters_1.SetMinSize( .005 )
# NETGEN_3D_Parameters_1.SetSecondOrder( 0 )
# NETGEN_3D_Parameters_1.SetOptimize( 1 )
# NETGEN_3D_Parameters_1.SetFineness( 4 )
# NETGEN_3D_Parameters_1.SetChordalError( -1 )
# NETGEN_3D_Parameters_1.SetChordalErrorEnabled( 0 )
# NETGEN_3D_Parameters_1.SetUseSurfaceCurvature( 1 )
# NETGEN_3D_Parameters_1.SetFuseEdges( 1 )
# NETGEN_3D_Parameters_1.SetQuadAllowed( 0 )
# NETGEN_3D_Parameters_1.SetCheckChartBoundary( 208 )
#
# NETGEN_1D_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=face_group0)
# NETGEN_2D_Parameters_1 = NETGEN_1D_2D_1.Parameters()
# NETGEN_2D_Parameters_1.SetMaxSize( 1)
# NETGEN_2D_Parameters_1.SetMinSize( .005 )
# NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
# NETGEN_2D_Parameters_1.SetOptimize( 1 )
# NETGEN_2D_Parameters_1.SetFineness( 4 )
# NETGEN_2D_Parameters_1.SetChordalError( -1 )
# NETGEN_2D_Parameters_1.SetChordalErrorEnabled( 0 )
# NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
# NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
# NETGEN_2D_Parameters_1.SetWorstElemMeasure( 0 )
# NETGEN_2D_Parameters_1.SetUseDelauney( 109 )
# NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
# NETGEN_2D_Parameters_1.SetCheckChartBoundary( 192 )
#
# NETGEN_1D_2D_2 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=face_group1)
# NETGEN_2D_Parameters_2 = NETGEN_1D_2D_2.Parameters()
# NETGEN_2D_Parameters_2.SetMaxSize( 1)
# NETGEN_2D_Parameters_2.SetMinSize( .005 )
# NETGEN_2D_Parameters_2.SetSecondOrder( 0 )
# NETGEN_2D_Parameters_2.SetOptimize( 1 )
# NETGEN_2D_Parameters_2.SetFineness( 4 )
# NETGEN_2D_Parameters_2.SetChordalError( -1 )
# NETGEN_2D_Parameters_2.SetChordalErrorEnabled( 0 )
# NETGEN_2D_Parameters_2.SetUseSurfaceCurvature( 1 )
# NETGEN_2D_Parameters_2.SetFuseEdges( 1 )
# NETGEN_2D_Parameters_2.SetWorstElemMeasure( 0 )
# NETGEN_2D_Parameters_2.SetUseDelauney( 109 )
# NETGEN_2D_Parameters_2.SetQuadAllowed( 0 )
# NETGEN_2D_Parameters_2.SetCheckChartBoundary( 192 )
#
# NETGEN_1D_2D_3 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D,geom=face_group2)
# NETGEN_2D_Parameters_3 = NETGEN_1D_2D_3.Parameters()
# NETGEN_2D_Parameters_3.SetMaxSize( 1 )
# NETGEN_2D_Parameters_3.SetMinSize( .005 )
# NETGEN_2D_Parameters_3.SetSecondOrder( 0 )
# NETGEN_2D_Parameters_3.SetOptimize( 1 )
# NETGEN_2D_Parameters_3.SetFineness( 4 )
# NETGEN_2D_Parameters_3.SetChordalError( -1 )
# NETGEN_2D_Parameters_3.SetChordalErrorEnabled( 0 )
# NETGEN_2D_Parameters_3.SetUseSurfaceCurvature( 1 )
# NETGEN_2D_Parameters_3.SetFuseEdges( 1 )
# NETGEN_2D_Parameters_3.SetWorstElemMeasure( 0 )
# NETGEN_2D_Parameters_3.SetUseDelauney( 109 )
# NETGEN_2D_Parameters_3.SetQuadAllowed( 0 )
# NETGEN_2D_Parameters_3.SetCheckChartBoundary( 192 )
#
# isDone = Mesh_1.Compute()
# Sub_mesh_1 = NETGEN_1D_2D_1.GetSubMesh()
# Sub_mesh_2 = NETGEN_1D_2D_2.GetSubMesh()
# Sub_mesh_3 = NETGEN_1D_2D_3.GetSubMesh()
#
# #
# airb1 = Mesh_1.GroupOnGeom(airb, 'airb', SMESH.VOLUME)
# airt1 = Mesh_1.GroupOnGeom(airt, 'airt', SMESH.VOLUME)
# OCA1 = Mesh_1.GroupOnGeom(OCA, 'OCA', SMESH.VOLUME)
# spine1 = Mesh_1.GroupOnGeom(spine, 'spine', SMESH.VOLUME)
# PVB1 = Mesh_1.GroupOnGeom(PVB, 'PVB', SMESH.VOLUME)
# sensor1 = Mesh_1.GroupOnGeom(sensor, 'sensor', SMESH.VOLUME)
#
# all_faces = [];
#
# for i in range(len(face0_names)):
#     all_faces.append(Mesh_1.GroupOnGeom(faces0[i], face0_names[i], SMESH.FACE))
# for i in range(len(face1_names)):
#     all_faces.append(Mesh_1.GroupOnGeom(faces1_new[i], face1_names[i], SMESH.FACE))
# for i in range(len(face2_names)):
#     all_faces.append(Mesh_1.GroupOnGeom(faces2_new[i], face2_names[i], SMESH.FACE))
#
# Mesh_1.ExportUNV('mesh.unv')
#
# Groups = Mesh_1.GetGroups()
#
# smesh.SetName(NETGEN_1D_2D_3D, 'NETGEN 1D-2D-3D')
# smesh.SetName(NETGEN_3D_Parameters_1, 'NETGEN 3D Parameters_1')
# smesh.SetName(NETGEN_1D_2D_1, 'NETGEN 1D-2D-1')
# smesh.SetName(NETGEN_2D_Parameters_1, 'NETGEN 3D Parameters_1')
# smesh.SetName(NETGEN_1D_2D_2, 'NETGEN 1D-2D-2')
# smesh.SetName(NETGEN_2D_Parameters_2, 'NETGEN 3D Parameters_2')
# smesh.SetName(NETGEN_1D_2D_3, 'NETGEN 1D-2D-3')
# smesh.SetName(NETGEN_2D_Parameters_3, 'NETGEN 3D Parameters_3')
#
# for i in range(6):
#     smesh.SetName(Groups[i],'Solid_' + str(i))
# for i in range(6,len(Groups)):
#     smesh.SetName(Groups[i],'Face_' + str(i-6))
#
# if salome.sg.hasDesktop():
#   salome.sg.updateObjBrowser()
