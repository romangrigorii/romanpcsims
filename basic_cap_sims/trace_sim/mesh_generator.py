#!/usr/bin/env python
# This file generates a series of FEM meshes which will then be used to compute capacitances for

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

distances = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3];
# distances = [.1]

cap_w = .03;
cap_l = 5;
N_lines = 8;

d_air = 4
d_glass = .4
d_spine = 0.7
d_OCA = 0.25
m = 0

for d in distances:

    geompy = geomBuilder.New()

    # generate geometrical entities

    O = geompy.MakeVertex(0, 0, 0)
    OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
    OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
    start = geompy.MakeFaceHW(cap_w, cap_l, 1)
    list_of_faces = [geompy.MakeTranslation(start, cap_w/2, cap_l/2, 0)];

    for i in range(N_lines-1):
        list_of_faces.append(geompy.MakeTranslation(list_of_faces[-1], cap_w + d, 0, 0))

    box_air0 = geompy.MakeBoxDXDYDZ(cap_w*N_lines + d*(N_lines-1), cap_l, d_air + d_glass + d_spine + d_OCA)
    box_air1 = geompy.MakeTranslation(box_air0,0,0, -(d_air + d_glass + d_spine + d_OCA)/2 - ((d_spine + d_OCA - d_glass)/2))

    box_glass0 = geompy.MakeBoxDXDYDZ(cap_w*N_lines + d*(N_lines-1),cap_l,d_glass)
    box_glass = geompy.MakeTranslation(box_glass0,0,0,0)

    box_spine0 = geompy.MakeBoxDXDYDZ(cap_w*N_lines + d*(N_lines-1),cap_l,d_spine)
    box_spine = geompy.MakeTranslation(box_spine0,0,0,-d_spine)

    box_OCA0 = geompy.MakeBoxDXDYDZ(cap_w*N_lines + d*(N_lines-1),cap_l,d_OCA)
    box_OCA = geompy.MakeTranslation(box_OCA0,0,0,-d_spine - d_OCA)

    box_air = geompy.MakeCutList(box_air1,[box_glass,box_spine,box_OCA],True)


    #identify each partition element
    list_of_entities = list_of_faces + [box_OCA,box_air,box_glass,box_spine]

    partition = geompy.MakePartition(list_of_entities, [], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
    [airb,OCA,spine,glass,airt] = geompy.ExtractShapes(partition, geompy.ShapeType["SOLID"], True)
    Faces = geompy.ExtractShapes(partition, geompy.ShapeType["FACE"], True)

    geompy.addToStudy( O, 'O' )
    geompy.addToStudy( OX, 'OX' )
    geompy.addToStudy( OY, 'OY' )
    geompy.addToStudy( OZ, 'OZ' )
    geompy.addToStudy( start, 'start' )
    for i in range(N_lines):
        geompy.addToStudy(list_of_faces[i], 'face' + str(i))

    geompy.addToStudy(airt,'airt');
    geompy.addToStudy(glass,'glass');
    geompy.addToStudy(spine,'spine');
    geompy.addToStudy(OCA,'OCA');
    geompy.addToStudy(airb,'airb');
    geompy.addToStudy(partition,'partition')
    geompy.addToStudyInFather(partition,airt,'airt');
    geompy.addToStudyInFather(partition,glass,'glass');
    geompy.addToStudyInFather(partition,spine,'spine');
    geompy.addToStudyInFather(partition,OCA,'OCA');
    geompy.addToStudyInFather(partition,airb,'airb');


    ### locating important Faces

    Vz = geompy.MakeVectorDXDYDZ(0, 0, 1)
    GND_vertex = geompy.MakeVertex(0, 0, -(d_spine+d_OCA))
    P_vertex = geompy.MakeVertex(0, 0, 0)
    airt_vertex = geompy.MakeVertex(0, 0, d_air/2 + d_glass);
    airb_vertex = geompy.MakeVertex(0, 0, -(d_air/2 + d_spine + d_OCA))

    face_group0 = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])
    face_group1 = geompy.CreateGroup(partition, geompy.ShapeType["FACE"])

    GND_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                        Vz, GND_vertex, GEOM.ST_ON)
    airt_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                        Vz, airt_vertex, GEOM.ST_ON)
    airb_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                        Vz, airb_vertex, GEOM.ST_ON)
    P_plane = geompy.GetShapesOnPlaneWithLocation(partition, geompy.ShapeType["FACE"],
                                                        Vz, P_vertex, GEOM.ST_ON)

    geompy.UnionList(face_group0,[airb_plane[0],GND_plane[0],airt_plane[0]])
    geompy.UnionList(face_group1,P_plane)

    geompy.addToStudyInFather(partition,face_group0,'face_group0')
    faces0 = geompy.ExtractShapes(face_group0, geompy.ShapeType["FACE"], True)
    face0_names = ['airb_face','GND_face','airt_face']
    for i in range(len(face0_names)):
        geompy.addToStudyInFather(face_group0,faces0[i],face0_names[i])

    geompy.addToStudyInFather(partition,face_group1,'face_group1')
    faces1 = geompy.ExtractShapes(face_group1, geompy.ShapeType["FACE"], True)
    face1_names = []
    face1_indexes = []
    print(len(faces1))
    for i in range(N_lines):
        face1_names.append('face' + str(i))
        face1_indexes.append(i*2)
        geompy.addToStudyInFather(face_group1,faces1[i*2],face1_names[i])

    ###
    ### SMESH component
    ###

    smesh = smeshBuilder.New()
    Mesh_1 = smesh.Mesh(partition)
    NETGEN_1D_2D_3D = Mesh_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    NETGEN_3D_Parameters_1 = NETGEN_1D_2D_3D.Parameters()
    NETGEN_3D_Parameters_1.SetMaxSize( .1 )
    NETGEN_3D_Parameters_1.SetMinSize( .001 )
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

    isDone = Mesh_1.Compute()
    Sub_mesh_1 = NETGEN_1D_2D_1.GetSubMesh()
    Sub_mesh_2 = NETGEN_1D_2D_2.GetSubMesh()

    airb1 = Mesh_1.GroupOnGeom(airb, 'airb', SMESH.VOLUME)
    airt1 = Mesh_1.GroupOnGeom(airt, 'airt', SMESH.VOLUME)
    OCA1 = Mesh_1.GroupOnGeom(OCA, 'OCA', SMESH.VOLUME)
    spine1 = Mesh_1.GroupOnGeom(spine, 'spine', SMESH.VOLUME)
    glass1 = Mesh_1.GroupOnGeom(glass, 'glass', SMESH.VOLUME)


    for i in range(len(face0_names)):
        Mesh_1.GroupOnGeom(faces0[i], face0_names[i], SMESH.FACE)
    for i in range(len(face1_names)):
        Mesh_1.GroupOnGeom(faces1[face1_indexes[i]], face1_names[i], SMESH.FACE)

    Mesh_1.ExportUNV(r'/media/sf_Shared_Folder/basic_cap_sims/trace_sim/UNV/Mesh_1_' + f"{int(d*1000):04}" + '.unv')

    print("mesh progress: " + str(round(100*m/len(distances))) + "%")
    m+=1

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
