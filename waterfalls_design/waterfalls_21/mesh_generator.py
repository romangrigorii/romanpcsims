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

maja = 3.2
mina = 2.2

Rx_d = 0.2

inum = 4
nlines = inum*2 + 1

step = 5.39

step_w = step
step_h = step/2

i_d = .385

cell_w_Hx = step_w - 2*delet - d_island
cell_w_Tx = step_w - 2*delet - d_island

d_air = 5
d_glass = 0.4
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

## CREATE FACES

# HX

gen_step_wx = geompy.MakeFaceHW(cell_w_Hx, 4*step_h, 1)

Hx0 = geompy.MakeTranslation(gen_step_wx, step_w/2, step_h*2, d_glass)
Hxlist = [geompy.MakeTranslation(gen_step_wx, step_w/2, step_h*2, d_glass)]
Hxlist.append(geompy.MakeTranslation(Hxlist[-1], step_w, 0, 0))

# IX

gen_i = geompy.MakeFaceHW(d_island, 4*step_h, 1)
Ixlist = [geompy.MakeTranslation(gen_i,step_w, step_h*2, d_glass)]
# for i in range(1):
#     Ixlist.append(geompy.MakeTranslation(Ixlist[-1],step_w, 0, 0))

# TXIX

TxIxlist = [geompy.MakeFaceHW(d_island, step_h - Rx_d - delet*2, 1)]
for i in range(4):
    TxIxlist.append(geompy.MakeTranslation(TxIxlist[-1],0,step_h,0))
TxIxlist = [geompy.MakeFuseList(TxIxlist,True,True)]
for h in range(3):
    TxIxlist.append(geompy.MakeTranslation(TxIxlist[0],step_w,0,0))

# IIX

base = geompy.MakeFaceHW(mina,d_island,1)
base = geompy.MakeTranslation(base,step_w/4,0,d_glass)
side = geompy.MakeFaceHW(d_island, step_h*4 - d_island*4 - delet*2, 1)
side = geompy.MakeTranslation(side,step_w/4,0,d_glass)
IIxlist = [side,base]
for i in range(inum):
    for p in [-1,1]:
        IIxlist.append(geompy.MakeTranslation(base,0,p*(i+1)*i_d,0))

for i in range(inum+1,inum*3 + 1,2):
    for p in [-1,1]:
        IIxlist.append(geompy.MakeTranslation(base,0,p*(i+1)*i_d,0))

IIxlist = geompy.MakeFuseList(IIxlist,True,True)
IIxlist = [geompy.MakeTranslation(IIxlist,0,-3.5*step_h,0)]
for i in range(2):
    IIxlist.append(geompy.MakeTranslation(IIxlist[-1],0,step_h*4,0))

IIxlist = [geompy.MakeFuseList(IIxlist,True,True)]
IIxlist.append(geompy.MakeTranslation(IIxlist[0],step_w/2,step_h*2,0));
IIxlist = [geompy.MakeFuseList(IIxlist,True,True)]
for h in range(1):
    if np.mod(h,2)==0:
        IIxlist.append(geompy.MakeTranslation(IIxlist[-1],step_w,step_h,0))
    else:
        IIxlist.append(geompy.MakeTranslation(IIxlist[-1],step_w,-step_h,0))

print(len(IIxlist))

### create new Hx

Hxlist_i = [];

base = geompy.MakeFaceHW(mina + delet*2,d_island + delet*2,1)
base = geompy.MakeTranslation(base,step_w/4,0,d_glass)
side = geompy.MakeFaceHW(d_island + delet*2, step_h*4 - d_island*4, 1)
side = geompy.MakeTranslation(side,step_w/4,0,d_glass)
IIxlist_d = [side,base]

for i in range(inum):
    for p in [-1,1]:
        IIxlist_d.append(geompy.MakeTranslation(base,0,p*(i+1)*i_d,0))

for i in range(inum+1,inum*3 + 1,2):
    for p in [-1,1]:
        IIxlist_d.append(geompy.MakeTranslation(base,0,p*(i+1)*i_d,0))

IIxlist_d = geompy.MakeFuseList(IIxlist_d,True,True)
IIxlist_d = [geompy.MakeTranslation(IIxlist_d,0,-3.5*step_h,0)]
for i in range(2):
    IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[-1],0,step_h*4,0))
IIxlist_d = [geompy.MakeFuseList(IIxlist_d,True,True)]
IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[0],step_w/2,step_h*2,0));
IIxlist_d = [geompy.MakeFuseList(IIxlist_d,True,True)]

for h in range(1):
    if np.mod(h,2)==0:
        IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[-1],step_w,step_h,0))
    else:
        IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[-1],step_w,-step_h,0))

for h in range(2):
    Hxlist_i.append(geompy.MakeCut(Hxlist[h], IIxlist_d[h], True))

Hxlist = Hxlist_i

# dummy islands

Hxlist_i = [];

base = geompy.MakeFaceHW(mina + delet*2,delet,1)
base1 = geompy.MakeTranslation(base,step_w/4,(delet+d_island)/2,d_glass)
base2 = geompy.MakeTranslation(base,step_w/4,-(delet+d_island)/2,d_glass)
base = geompy.MakeFuseList([base1,base2],True,True)
IIxlist_d = [base]

for i in range(inum*4 + 1):
    for p in [-1,1]:
        IIxlist_d.append(geompy.MakeTranslation(base,0,p*(i+1)*i_d,0))

IIxlist_d = geompy.MakeFuseList(IIxlist_d,True,True)
IIxlist_d = [geompy.MakeTranslation(IIxlist_d,0,-3.5*step_h,0)]
for i in range(2):
    IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[-1],0,step_h*4,0))
IIxlist_d = [geompy.MakeFuseList(IIxlist_d,True,True)]
IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[0],step_w/2,step_h*2,0));
IIxlist_d = [geompy.MakeFuseList(IIxlist_d,True,True)]

for h in range(1):
    if np.mod(h,2)==0:
        IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[-1],step_w,step_h,0))
    else:
        IIxlist_d.append(geompy.MakeTranslation(IIxlist_d[-1],step_w,-step_h,0))

for h in range(2):
    Hxlist_i.append(geompy.MakeCut(Hxlist[h], IIxlist_d[h], True))

Hxlist = Hxlist_i

## creating Rx

gen_rx_f = geompy.MakeFaceHW(mina, maja, 1)
gen_rx_t = geompy.MakeFaceHW(step_w, Rx_d, 1)
gen_rx_ft = geompy.MakeFuseList([gen_rx_f, gen_rx_t], True, True)


Rx0_0 = geompy.MakeTranslation(gen_rx_ft,step_w*1/4, step_h/2, 0)
Rx0_1 = geompy.MakeTranslation(gen_rx_t,step_w*1/4, step_h*3/2, 0)
Rx0_2 = geompy.MakeTranslation(gen_rx_t,step_w*1/4, step_h*5/2, 0)
Rx0_3 = geompy.MakeTranslation(gen_rx_t,step_w*1/4, step_h*7/2, 0)

Rx1_0 = geompy.MakeTranslation(gen_rx_t,step_w*3/4, step_h/2, 0)
Rx1_1 = geompy.MakeTranslation(gen_rx_t,step_w*3/4, step_h*3/2, 0)
Rx1_2 = geompy.MakeTranslation(gen_rx_ft,step_w*3/4, step_h*5/2, 0)
Rx1_3 = geompy.MakeTranslation(gen_rx_t,step_w*3/4, step_h*7/2, 0)

Rx2_0 = geompy.MakeTranslation(gen_rx_t,step_w*5/4, step_h/2, 0)
Rx2_1 = geompy.MakeTranslation(gen_rx_ft,step_w*5/4, step_h*3/2, 0)
Rx2_2 = geompy.MakeTranslation(gen_rx_t,step_w*5/4, step_h*5/2, 0)
Rx2_3 = geompy.MakeTranslation(gen_rx_t,step_w*5/4, step_h*7/2, 0)

Rx3_0 = geompy.MakeTranslation(gen_rx_t,step_w*7/4, step_h/2, 0)
Rx3_1 = geompy.MakeTranslation(gen_rx_t,step_w*7/4, step_h*3/2, 0)
Rx3_2 = geompy.MakeTranslation(gen_rx_t,step_w*7/4, step_h*5/2, 0)
Rx3_3 = geompy.MakeTranslation(gen_rx_ft,step_w*7/4, step_h*7/2, 0)

Rx0 = geompy.MakeFuseList([Rx0_0, Rx1_0, Rx2_0, Rx3_0], True, True)
Rx1 = geompy.MakeFuseList([Rx0_1, Rx1_1, Rx2_1, Rx3_1], True, True)
Rx2 = geompy.MakeFuseList([Rx0_2, Rx1_2, Rx2_2, Rx3_2], True, True)
Rx3 = geompy.MakeFuseList([Rx0_3, Rx1_3, Rx2_3, Rx3_3], True, True)

Rxlist = [Rx0, Rx1, Rx2, Rx3]

# TX

gen_rx_deletion_f = geompy.MakeFaceHW(mina + delet*2,maja + delet*2, 1)
gen_rx_deletion_fL = geompy.MakeTranslation(gen_rx_deletion_f,-step_w/4,0,0)
gen_rx_deletion_fR = geompy.MakeTranslation(gen_rx_deletion_f,step_w/4,0,0)
gen_rx_deletion_t = geompy.MakeFaceHW(step_w, Rx_d + delet*2, 1)
gen_rx_deletion_ftL = geompy.MakeFuseList([gen_rx_deletion_fL,gen_rx_deletion_t])
gen_rx_deletion_ftR = geompy.MakeFuseList([gen_rx_deletion_fR,gen_rx_deletion_t])

gen_cell_Tx = geompy.MakeFaceHW(cell_w_Tx, step_h, 1)

gen_tx_0L = geompy.MakeCutList(gen_cell_Tx, [gen_rx_deletion_ftL], True)
gen_tx_0R = geompy.MakeCutList(gen_cell_Tx, [gen_rx_deletion_ftR], True)

gen_tx_1L = geompy.MakeCutList(gen_cell_Tx, [gen_rx_deletion_t,geompy.MakeTranslation(gen_rx_deletion_fL,0,-step_h,0),geompy.MakeTranslation(gen_rx_deletion_fR,0,step_h,0)], True)
gen_tx_1R = geompy.MakeCutList(gen_cell_Tx, [gen_rx_deletion_t,geompy.MakeTranslation(gen_rx_deletion_fR,0,-step_h,0),geompy.MakeTranslation(gen_rx_deletion_fL,0,step_h,0)], True)

Tx0_0 = geompy.MakeTranslation(gen_tx_0L,step_w*1/2, step_h/2, 0)
Tx0_1 = geompy.MakeTranslation(gen_tx_1L,step_w*1/2, step_h*3/2, 0)
Tx0_2 = geompy.MakeTranslation(gen_tx_0R,step_w*1/2, step_h*5/2, 0)
Tx0_3 = geompy.MakeTranslation(gen_tx_1R,step_w*1/2, step_h*7/2, 0)

Tx1_0 = geompy.MakeTranslation(gen_tx_1R,step_w*3/2, step_h/2, 0)
Tx1_1 = geompy.MakeTranslation(gen_tx_0L,step_w*3/2, step_h*3/2, 0)
Tx1_2 = geompy.MakeTranslation(gen_tx_1L,step_w*3/2, step_h*5/2, 0)
Tx1_3 = geompy.MakeTranslation(gen_tx_0R,step_w*3/2, step_h*7/2, 0)

Tx0 = geompy.MakeFuseList([Tx0_0, Tx0_1, Tx0_2, Tx0_3], True, True)
Tx1 = geompy.MakeFuseList([Tx1_0, Tx1_1, Tx1_2, Tx1_3], True, True)

Txlist = [Tx0, Tx1]

### create solids

box_air_surr = geompy.MakeBoxDXDYDZ(step_w*2,step_h*4,d_air)
box_air1 = geompy.MakeTranslation(box_air_surr, 0, 0, d_glass)
box_air2 = geompy.MakeTranslation(box_air_surr, 0, 0, -(d_PVB+d_spine+d_OCA+d_air))

box_glass0 = geompy.MakeBoxDXDYDZ(step_w*2,step_h*4,d_glass)
box_glass = geompy.MakeTranslation(box_glass0,0,0,0)

box_PVB0 = geompy.MakeBoxDXDYDZ(step_w*2,step_h*4,d_PVB)
box_PVB = geompy.MakeTranslation(box_PVB0,0,0,-d_PVB)

box_spine0 = geompy.MakeBoxDXDYDZ(step_w*2,step_h*4,d_spine)
box_spine = geompy.MakeTranslation(box_spine0,0,0,-d_PVB - d_spine)

box_OCA0 = geompy.MakeBoxDXDYDZ(step_w*2,step_h*4,d_OCA)
box_OCA = geompy.MakeTranslation(box_OCA0,0,0,-d_PVB - d_spine - d_OCA)

# # making the partitions

list_of_bodies = Rxlist + Txlist + Hxlist + Ixlist + IIxlist + TxIxlist + [box_OCA,box_PVB,box_air1,box_air2,box_glass,box_spine]
partition = geompy.MakePartition(list_of_bodies, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
[airb,OCA,spine,PVB,glass,airt] = geompy.ExtractShapes(partition, geompy.ShapeType["SOLID"], True)
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
for i in range(len(TxIxlist)):
    geompy.addToStudy( TxIxlist[i], 'TxIx' + str(i))

geompy.addToStudy( box_air_surr, 'air_s' )
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
airt_vertex = geompy.MakeVertex(0, 0, d_air + d_glass);
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

#### identifying the faces and grouping, Hx layer

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
# face1_indexes_Hx = [5, 15,19, 29];
# face1_names_Ix = [];
# face1_indexes_Ix = [11,23];
# face1_names_IIx = [];
# face1_indexes_IIx = [2,4,6,8,  13,17,21, 26,28,30,32];
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
    if (np.abs(xmi-xma)< d_island*1.05) and (np.abs(xmi-xma) > d_island*.95):
            face1_indexes_Ix.append(i)
            face1_names_Ix.append('Ix' + str(ixi))
            ixi+=1
    if (np.abs(xmi - xma) < mina + delet) and (np.abs(xmi - xma) > mina - delet):
        if (np.abs(ymi - yma) >  delet*2):
            face1_indexes_IIx.append(i)
            face1_names_IIx.append('IIx' + str(iixi))
            iixi+=1
    if (np.abs(xmi - xma) <  cell_w_Hx*1.05) and (np.abs(xmi - xma) > cell_w_Hx*.95):
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

#### identifying the faces and grouping, Tx layer

geompy.addToStudyInFather(partition,face_group2,'face_group2')
faces2 = geompy.ExtractShapes(face_group2, geompy.ShapeType["FACE"], True)

# face2_names_TxIx = [];
# face2_indexes_TxIx = [];
# face2_names_Tx = [];
# face2_indexes_Tx = [];
# face2_names_Rx = [];
# face2_indexes_Rx = [];
# for i in range(len(faces2)):
#     face2_names_Tx.append(str(i))
#     face2_indexes_Tx.append(i)

face2_names_TxIx = [];
face2_indexes_TxIx = [20,21,25,28,30]
face2_names_Tx = [];
face2_indexes_Tx = [1,4,8,9,12,13,14,15,16,17,23,26, 24,27,34,35,36,37,38,39,42,43,44,45]
face2_names_Rx = [];
face2_indexes_Rx = [11,33,18,40];
for i in range(len(face2_indexes_Tx)):
    face2_names_Tx.append("Tx" + str(i))
for i in range(len(face2_indexes_Rx)):
    face2_names_Rx.append("Rx" + str(i))
for i in range(len(face2_indexes_TxIx)):
    face2_names_TxIx.append("TxIx" + str(i))

# face2_names_TxIx = [];
# face2_indexes_TxIx = [];
# face2_names_Tx = [];
# face2_indexes_Tx = [];
# face2_names_Rx = [];
# face2_indexes_Rx = [];
# txi = 0;
# rxi = 0;
# tixi = 0;
# for i in range(len(faces1)):
#     [xmi,xma,ymi,yma,zmi,zma] = geompy.BoundingBox(faces1[i])
#     if (np.abs(xmi-xma)< d_island*1.05) and (np.abs(xmi-xma) > d_island*.95):
#         face2_indexes_TxIx.append(i)
#         face2_names_TxIx.append('TxIx' + str(tixi))
#         tixi+=1
#     if (np.abs(xmi - xma) < cell_w_Hx*1.05) and (np.abs(xmi - xma) > cell_w_Hx/5):
#         face2_indexes_Tx.append(i)
#         face2_names_Tx.append('Tx' + str(txi))
#         txi+=1
#     if (np.abs(xmi - xma) <  cell_w_Hx*3*1.05) and (np.abs(xmi - xma) > cell_w_Hx*3*.95):
#         face2_indexes_Rx.append(i)
#         face2_names_Rx.append('Rx' + str(rxi))
#         rxi+=1

face2_names = face2_names_Tx + face2_names_Rx + face2_names_TxIx
face2_indexes = face2_indexes_Tx + face2_indexes_Rx + face2_indexes_TxIx
faces2_new = [];
for i in face2_indexes:
    faces2_new.append(faces2[i])

for i in range(len(face2_names)):
    geompy.addToStudyInFather(face_group2,faces2_new[i],face2_names[i])

#### create mesh groups
#
# smesh = smeshBuilder.New()
# Mesh_1 = smesh.Mesh(partition)
# NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
# NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
# NETGEN_3D_Parameters_1.SetMaxSize( 1 )
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
# glass1 = Mesh_1.GroupOnGeom(glass, 'glass', SMESH.VOLUME)
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
