#!/usr/bin/env python3
# -*- coding: utf-8 -*-



from scipy.spatial import Voronoi, voronoi_plot_2d
from Graph import Graph, addLimitPoints
from ForceDirectedLayout import ForceDirectedLayout
import numpy as np

#############################################################

points = np.loadtxt('points2Dtxt.csv')
edgesCoords = np.load('edgesCoords.npy')
braidCrossings = np.load('braidsInfo.npy')
points = addLimitPoints(points)
vor = Voronoi(points)
#voronoi_plot_2d(vor)
vertices = vor.vertices

def dist(v1,v2):
    return np.sqrt((v1[0]-v2[0])**2+(v1[1]-v2[1])**2)

def findNearest(v,vertices):
    distance = -1
    candidate = -1
    for i in range(len(vertices)):
        d = dist(v,vertices[i])
        if distance == -1 or distance > d:
            distance = d
            candidate = i
    return candidate

braidEdges =[]
for e in edgesCoords:
    v0 = findNearest(e[0],vertices)
    v1 = findNearest(e[1],vertices)
    e_new = [v0,v1] if v0<v1 else [v1,v0]
    braidEdges.append(e_new)

braidInfo = [braidEdges, braidCrossings]

g = Graph(vor,braidInfo)
g.plot_graph()
g.plot = True
# Force algorithm
maxIter = 107
# maxIter = 2
tol = 0.2
beta = 1000*g.calculate_scale()  # node_node_attraction
delta = 0.001  # node_node_repulsion
gamma = 0.01 # node_edge_repulsion
f = ForceDirectedLayout(g, beta, delta, gamma, maxIter, opt=True)
f.run()
