#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 22:42:06 2019

@author: calquezar
"""
import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

#########################################################
# Necesidad de añadir cuatro puntos para cerrar el diagrama de voronoi
points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2],\
                   [-5,5],[5,5],[5,-5],[-5,-5]])
vor = Voronoi(points)
#voronoi_plot_2d(vor)
#plt.show()
#########################################################
# Input para el algoritmo ImPrEd
vertices = vor.vertices
edges = [ x for x in vor.ridge_vertices if -1 not in x ]
regions = [ x for x in vor.regions if -1 not in x ][1:] # hay que asegurarse que el primer elementos es [] siempre
#########################################################
# Creamos la clase Sv

class Sv:
  def __init__(self,v,e,c):
    self.vertices = v
    self.edges = e
    self.connectedTo = c
#########################################################
# Calculamos la estructura Sv
# Todos los vertices y aristas de Sv para todo v

allSv =[]
for v in range(len(vertices)):
  #Faces
  vFaces = []
  for r in regions:
    if v in r:
      vFaces += r
#  vFaces = list(filter(lambda x: x != v, vFaces))
  vFaces = list(set(vFaces)) # remove redundancies
  vFaces.remove(v)
  vFaces.sort()
  print(vFaces)
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
#########################################################
  
  
