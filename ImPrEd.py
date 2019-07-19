#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 22:42:06 2019

@author: calquezar
"""
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import LineString, Point
import matplotlib.pyplot as plt
import math
from  QuadTree import QPoint, QuadTree
#########################################################
# Sv class
class Sv:
  def __init__(self,v,e,c):
    self.vertices = v
    self.edges = e
    self.connectedTo = c
#########################################################
def preprocessVoronoiStruct(vor):
  vertices = vor.vertices
  edges = [ x for x in vor.ridge_vertices if -1 not in x ]
  regions = [ x for x in vor.regions if -1 not in x ][1:] # hay que asegurarse que el primer elementos es [] siempre
  return vertices,edges,regions
#########################################################
#Forces
def dist(v1,v2):
  return math.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2)
# Fuerza de repulsion entre vertices
def fRep(v1,v2,delta=1):
  d = dist(v1,v2)
  fx = ((delta/d)**2)*(v1[0]-v2[0])
  fy = ((delta/d)**2)*(v1[1]-v2[1])
  return [fx,fy]
# Fuerza de atracción entre vertices conectados
def fattract(v1,v2,delta=1):
  d = dist(v1,v2)
  fx = (d/delta)*(v2[0]-v1[0])
  fy = (d/delta)*(v2[1]-v1[1])
  return [fx,fy]
# Fuerza de repulsión entre vertice y arista. Ve es la proyección de v sobre e
def fvEdge(v,ve,gamma=1):
  d = dist(v,ve)
  fx=0
  fy=0
  if(d<gamma):
   fx = (((gamma-d)**2)/d)*(v[0]-ve[0])
   fy = (((gamma-d)**2)/d)*(v[1]-ve[1])
  return [fx,fy]
#########################################################
# Calculamos la estructura Sv
# Todos los vertices y aristas de Sv para todo v
def preprocessing(vertices,edges,regions):
  allSv =[]
  for v in range(len(vertices)):
    #Faces
    vFaces = []
    for r in regions:
      if v in r:
        vFaces += r
  #  vFaces = list(filter(lambda x: x != v, vFaces))
    vFaces = list(set(vFaces)) # remove redundancies
    vFaces.remove(v) # remove myself from the list
    vFaces.sort() # sort ascending
    # Connections
    connectedTo = []
    for v2 in vFaces:
      e = [v,v2] if v<v2 else [v2,v]
      if(e in edges):
        connectedTo += [v2]
    connectedTo.sort()
    # Surrounding Edges
    eFaces = []
    for v1 in range(len(vFaces)-1):
      for v2 in range(v1+1,len(vFaces)):
        if [vFaces[v1],vFaces[v2]] in edges:#asumo las aristas se dan con los vertices ordenados de menor a mayor
          eFaces.append([vFaces[v1],vFaces[v2]])
    allSv.append(Sv(vFaces,eFaces,connectedTo))
  return allSv
  
#########################################################
def pointEdgeProjection(v,e): #TODO REVISAR
  line = LineString([e[0],e[1]])
  p = Point(v[0],v[1])
  nearestPoint = list(line.interpolate(line.project(p)).coords[0]) #using shapely to find the nearest point to line
  if (dist(nearestPoint,e[0])<dist(e[0],e[1])): #check if the point is inside the segment or not
    ve = nearestPoint
  else: # get the nearest end
    ve = e[0] if (dist(nearestPoint,e[0]) < dist(nearestPoint,e[1])) else e[1]
  return ve
#########################################################
# 1) Forces calculation for each vertex
def forcesCalculation(vertices,qTree,allSv,delta,gamma):
  forces = []
  for i in range(len(vertices)):
#  for i in range(1):
    v = vertices[i]
    sv = allSv[i]
    fa = [0,0]
    fr = [0,0]
    fe = [0,0]
    # Atracción entre nodos conectados
    for vIdx in sv.connectedTo:
      vc = vertices[vIdx]
      f = fattract(v,vc)
      fa[0]+=f[0]
#      fa[1]+=f[1]
    # Repulsion entre nodos
    neighbours = qTree.findPoints(QPoint(v[0],v[1]),delta)
    for n in neighbours:
      f = fRep(v,n.toArray())
      fr[0]+=f[0]
      fr[1]+=f[1]
    # Repulsión entre nodos y aristas
    for edge in sv.edges:
      edgeCoords = [vertices[edge[0]], vertices[edge[1]]] 
      ve = pointEdgeProjection(v,edgeCoords)
      f = fvEdge(v,ve,gamma)
      fe[0]+=f[0]
      fe[1]+=f[1]
    # Cálculo de la resultante de todas fuerzas
    fTotal = [0,0]
    fTotal[0] = fa[0]+fr[0]+fe[0]
    fTotal[1] = fa[1]+fr[1]+fe[1]
    forces.append(fTotal)
  return forces
#########################################################
# 2) Cálculo de Mv
def calculateMaxDisplacements(v,e):
  w1 = e[0]
  w2 = e[1]
  w = [vertices[w1],vertices[w2]]
  ve = pointEdgeProjection(v,w)
  cv = dist(v,ve)
  angleCv = math.atan2(ve[0]-v[0],ve[1]-v[1])+math.pi #angle between [0,2*pi]
  sectorCv = angleCv//(math.pi/4)
  distances = []
  for j in range(8):
    diff = (sectorCv-j)%8
    if diff==0: #same sector
      sigma=1
    elif (diff==1 or diff==2):
      sigma= 1/math.cos(angleCv-(j+1)*math.pi/4)
    elif (diff==6 or diff==7):
      sigma= 1/math.cos(angleCv-j*math.pi/4)
    else:
      sigma=2 # own decision
    distances.append(cv*sigma)
  return distances

def updateMvs(mv,mvNew):
  if (mv==[]):
    return mvNew
  else:
    for i in range(len(mv)):
      mv[i]=min(mv[i],mvNew[i])
    return mv

def calculateMvs(vertices,allSv):
  Mvs = []
  for i in range(len(vertices)):
    v = vertices[i]
    sv = allSv[i]
    mv = []
    for e in sv.edges:
      maxDisp = calculateMaxDisplacements(v,e)
      mv = updateMvs(mv,maxDisp)
    Mvs.append(mv)
  return Mvs
  
#########################################################
# 3) Desplazamiento de los vértices en base al min(F,Mv)
def moveNodes(vertices,forces,Mvs):
  for i in range(len(vertices)):
    angleF = math.atan2(forces[i][1],forces[i][0])+math.pi #angle between [0,2*pi]
    sector = int(angleF//(math.pi/4))%8
    forceAmplitude = math.sqrt(forces[i][0]**2+forces[i][1]**2)
    mv = Mvs[i][sector]
    maxDisp = min(forceAmplitude,mv)
    dispX = maxDisp*math.cos(angleF)
    dispY = maxDisp*math.sin(angleF)
    vertices[i][0]+=dispX
    vertices[i][1]+=dispY
  return vertices
  
def plotGraph(vertices,edges):
  plt.clf()
  for e in edges:
    v0 = vertices[e[0]]
    v1 = vertices[e[1]]
    x = [v0[0],v1[0]]
    y = [v0[1],v1[1]]
    plt.plot(x,y,color='k')
  axes = plt.gca()
  vmin = np.amin(vertices)
  vmax = np.amax(vertices)
  axes.set_xlim([vmin,vmax])
  axes.set_ylim([vmin,vmax])
  plt.pause(1)
#########################################################
# Datos de prueba
# Necesidad de añadir cuatro puntos para cerrar el diagrama de voronoi
#points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
#                   [2, 1], [2, 2], [-5,5],[5,5],[5,-5],[-5,-5]])
#vor = Voronoi(points)
#voronoi_plot_2d(vor)
#plt.show()
#########################################################
points = np.array([[1,1], [3, 1], [2, 2], [2,0],[1.7,1],[2.3,1 ],[2,1.7]])
vor = Voronoi(points)
vor.vertices[0][0]=1.7
vor.vertices[0][1]=1.3
#voronoi_plot_2d(vor)
#plt.show()
######################################################### #TODO
# Main algorithm
delta = 0.5
gamma = 1
vertices,edges,regions = preprocessVoronoiStruct(vor)
allSv = preprocessing(vertices,edges,regions)
maxIter = 4
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
for it in range(maxIter):
  qTree = QuadTree(QPoint.arrayToList(vertices))
#  qTree.plot()
  #Step 1
  forces = forcesCalculation(vertices,qTree,allSv,delta,gamma)
#  #Step 2
  Mvs = calculateMvs(vertices,allSv)
#  #Step 3
  moveNodes(vertices,forces,Mvs)
  #refresh plot
  plotGraph(vertices,edges)

plt.show()
  
  
#  plt.scatter(vertices[:,0], vertices[:,1],s=10)
#  plt.pause(0.1)


#Draw results
#vor.vertices = vertices
#voronoi_plot_2d(vor)
#plt.show()
#  

