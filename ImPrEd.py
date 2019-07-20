#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 22:42:06 2019

@author: calquezar
"""
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import LineString, Point, Polygon
import matplotlib.pyplot as plt
import math
from  QuadTree import QPoint, QuadTree
from statistics import median
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
  if(d>0):
    fx = ((delta/d)**2)*(v1[0]-v2[0])
    fy = ((delta/d)**2)*(v1[1]-v2[1])
    return [fx,fy]
  else:
    epsilon = np.finfo(np.float32).eps
    return [1/epsilon,1/epsilon]
# Fuerza de atracción entre vertices conectados
def fattract(v1,v2,delta=1):
  d = dist(v1,v2)
  fx = (d/delta)*(v2[0]-v1[0])
  fy = (d/delta)*(v2[1]-v1[1])
  return [fx,fy]
# Fuerza de repulsión entre vertice y arista. Ve es la proyección de v sobre e
def fvEdge(v,ve,gamma=1):
  d = dist(v,ve)
  if(d>0):
    fx=0
    fy=0
    if(d<gamma):
     fx = (((gamma-d)**2)/d)*(v[0]-ve[0])
     fy = (((gamma-d)**2)/d)*(v[1]-ve[1])
    return [fx,fy]
  else:
    epsilon = np.finfo(np.float32).eps
    return [1/epsilon,1/epsilon]
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
def forcesCalculation(vertices,region,qTree,allSv,delta,gamma):
  forces = []
#  for i in range(len(vertices)):
  for i in region:
#  for i in range(1):
    v = vertices[i]
    sv = allSv[i]
    fa = [0,0]
    fr = [0,0]
    fe = [0,0]
    # Atracción entre nodos conectados
    for vIdx in sv.connectedTo:
      vc = vertices[vIdx]
      f = fattract(v,vc,delta)
      fa[0]+=f[0]
      fa[1]+=f[1]
#    print(fa)
#     Repulsion entre nodos
    neighbours = qTree.findPoints(QPoint(v[0],v[1]),delta)
    for n in neighbours:
#      print(n.toArray())
      f = fRep(v,n.toArray(),delta)
      fr[0]+=f[0]
      fr[1]+=f[1]
#    print(fr)
    # Repulsión entre nodos y aristas
    for edge in sv.edges:
      edgeCoords = [vertices[edge[0]], vertices[edge[1]]] 
      ve = pointEdgeProjection(v,edgeCoords)
      f = fvEdge(v,ve,gamma)
      fe[0]+=f[0]
      fe[1]+=f[1]
    # Cálculo de la resultante de todas fuerzas
#    print(fa,fr,fe)
    fTotal = [0,0]
    fTotal[0] = 2*fa[0]+fr[0]+fe[0]
    fTotal[1] = 2*fa[1]+fr[1]+fe[1]
    forces.append(fTotal)
#    print(fTotal)
  return forces
#########################################################
# 2) Cálculo de Mv
def calculateMaxDisplacements(v,e):
  w1 = e[0]
  w2 = e[1]
  w = [vertices[w1],vertices[w2]]
  ve = pointEdgeProjection(v,w)
  cv = dist(v,ve)
  angleCv = (math.atan2(ve[1]-v[1],ve[0]-v[0])+2*math.pi)%(2*math.pi) #angle between [0,2*pi]
  sectorCv = int(angleCv//(math.pi/4))
  distances = []
  for j in range(8):
    diff = (sectorCv-j)%8
    if diff==0: #same sector
      sigma=1
    elif (diff==1 or diff==2):
      denominator=math.cos(angleCv-(j+1)*math.pi/4)
      if(denominator!=0):
        sigma= 1/denominator
      else:
        epsilon = np.finfo(np.float32).eps
        sigma=1/epsilon
    elif (diff==6 or diff==7):
      denominator=math.cos(angleCv-j*math.pi/4)
      if(denominator!=0):
        sigma= 1/denominator
      else:
        epsilon = np.finfo(np.float32).eps
        sigma=1/epsilon
    else:
      sigma=10 # own decision
    distances.append(cv*sigma/2)
  return distances

def updateMvs(mv,mvNew):
  if (mv==[]):
    return mvNew
  else:
    for i in range(len(mv)):
      mv[i]= mv[i] if abs(mv[i])<abs(mvNew[i]) else mvNew[i]
    return mv

def calculateMvs(vertices,region,allSv):
  Mvs = []
#  for i in range(len(vertices)):
  for i in region:
    v = vertices[i]
    sv = allSv[i]
    mv = []
    for e in sv.edges:
      maxDisp = calculateMaxDisplacements(v,e)
      mv = updateMvs(mv,maxDisp)
#    if(i==4):
#      print(mv)
    Mvs.append(mv)
  return Mvs
  
#########################################################
# 3) Desplazamiento de los vértices en base al min(F,Mv)
def moveNodes(vertices,region,forces,Mvs,verticesBoundary):
#  for i in range(len(vertices)):
  for i in range(len(region)):
    v = region[i]
    angleF = (math.atan2(forces[i][1],forces[i][0])+2*math.pi)%(2*math.pi) #angle between [0,2*pi]
    sector = int(angleF//(math.pi/4))%8
#    print(forces[i],angleF*180/math.pi,sector)
    forceAmplitude = math.sqrt(forces[i][0]**2+forces[i][1]**2)
    mv = Mvs[i][sector]
    maxShift = forceAmplitude if abs(forceAmplitude)<abs(mv) else mv
#    if(i==4):
#      print(forceAmplitude,mv,maxDisp)
    shiftX = maxShift*math.cos(angleF)
    shiftY = maxShift*math.sin(angleF)
#    print(shiftX,shiftY)
    vertices[v][0]+=shiftX
    vertices[v][1]+=shiftY
    plotGraph(vertices,edges)

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
  axes.set_xlim([vmin-0.1,vmax+0.1])
  axes.set_ylim([vmin-0.1,vmax+0.1])
  plt.pause(0.1)
#########################################################
  
  
def area(region):
  points =[vertices[v] for v in region]
  pol = Polygon(points)
  return pol.area
  
#########################################################
  
def isInRegion(edge,region):
  isIn = True
  for v in edge:
    if v not in region:
      isIn=False
      break
  return isIn

def getBoundary(edges,regions):
  
  boundary = []
  for e in edges:
    count = 0
    for r in regions:
      if isInRegion(e,r):
        count+=1
      if count>1:
        break
    if count<2:
      boundary.append(e)
  return boundary
  
#########################################################
# Datos de prueba
# Necesidad de añadir cuatro puntos para cerrar el diagrama de voronoi
points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
                   [2, 1], [2, 2], [-5,5],[5,5],[5,-5],[-5,-5]])
vor = Voronoi(points)
vor.vertices[13][0]=0.0000000005
vor.vertices[13][1]=0.0000000005
vor.vertices[14][0]=0.0000000015
vor.vertices[14][1]=0.0000000005
vor.vertices[15][0]=0.0000000015
vor.vertices[15][1]=0.0000000015
vor.vertices[12][0]=0.0000000005
vor.vertices[12][1]=0.0000000015

#voronoi_plot_2d(vor)
#plt.show()
#########################################################
#points = np.array([[-10,0],[10,0],[0,-10],[0,10],[-0.2,-0.5],[0, 0], [0.0001, 0.0001], [0, 0.0001],[0.0001, 0],[0.00005, 0.00005],[0, 0.0002],[0.0001, 0.0002]])
#vor = Voronoi(points)
#vor.vertices[2][0]=-10000000
#vor.vertices [2][1]=-10000000
#########################################################
#points = np.array([[1,1], [3, 1], [2, 2], [2,0],[1.7,1],[2.3,1 ],[2,1.7]])
#vor = Voronoi(points)
#vor.vertices[0][0]=1.7
#vor.vertices[0][1]=1.3
#voronoi_plot_2d(vor)
#plt.show()
#########################################################
#vertices = [[0,0],[10,0],[0.1,0.1],[0,10],[10,10]]
#edges = [[0,1],[0,2],[0,3],[1,2],[1,4],[2,3],[2,4],[3,4]]
#regions = [[0,1,2],[0,2,3],[1,2,4],[2,3,4]]
######################################################### #TODO
# Main algorithm
delta = 1
gamma = 1
vertices,edges,regions = preprocessVoronoiStruct(vor)
boundary = getBoundary(edges,regions)

verticesBoundary = []
for e in boundary:
  verticesBoundary+=e
verticesBoundary = set(verticesBoundary)
    
allSv = preprocessing(vertices,edges,regions)
maxIter = 20
figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
#plotGraph(vertices,edges)

for it in range(maxIter):
  print(it)
  qTree = QuadTree(QPoint.arrayToList(vertices))
#  qTree.plot()
  #####################################################
  # region selection in function of relative area respect to the median
  # of the graph
  regions.sort(key=area) # sort regions by area in ascending order
  areas = [area(r) for r in regions]
  mAreas = median(areas)
  minArea = areas[0]
  maxArea = areas[-1]
  region = regions[0] if abs(math.log(mAreas/minArea))>abs(math.log(mAreas/maxArea)) else regions[-1]
  #####################################################
  #Step 1
  forces = forcesCalculation(vertices,region,qTree,allSv,delta,gamma)
#  #Step 2
  Mvs = calculateMvs(vertices,region,allSv)
#  #Step 3
  moveNodes(vertices,region,forces,Mvs,verticesBoundary)
  #refresh plot
#plotGraph(vertices,edges)
plt.show()
  
  
#  plt.scatter(vertices[:,0], vertices[:,1],s=10)
#  plt.pause(0.1)


#Draw results
#vor.vertices = vertices
#voronoi_plot_2d(vor)
#plt.show()
#  

