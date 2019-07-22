#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Graph data structure
...description...

REFERENCES:
- :wikipedia:`Graph theory`

AUTHORS:
- Carlos Alquézar Baeta
"""

# ****************************************************************************
#       Copyright (C) 2019  Carlos Alquézar Baeta
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import numpy as np
from matplotlib import pyplot as plt
from shapely.geometry import Polygon

class SurroundingInfo:
    def __init__(self, v, e, c):
        self.vertices = v # vertices in the boundary of the surrounding area
        self.edges = e # edges forming the boundary of the surrounding area
        self.connected_to = c # edges that connect the vertex to others in the boundary
#########################################################

class Graph:

    def __init__(self, vor):
        r"""
            Constructor
        """
        self.vertices = vor.vertices
        self.edges = [ x for x in vor.ridge_vertices if -1 not in x ]
        self.regions = [ x for x in vor.regions if -1 not in x ][1:] # hay que asegurarse que el primer elementos es [] siempre
        self.plot = False

    def area(self, region):
        r"""
            Return the area of the region
        """
        points = [self.vertices[v] for v in region]
        pol = Polygon(points)
        return pol.area

    def is_in_region(edge, region):
        r"""
            Return true if the edge is in the boundary of the region
        """
        isIn = True
        for v in edge:
            if v not in region:
                isIn = False
                break
        return isIn

    def get_boundary_edges(self):
        r"""
            Return
        """
        boundary = []
        for e in self.edges:
            count = 0
            for r in self.regions:
                if Graph.is_in_region(e, r):
                    count += 1
                if count > 1:
                    break
            if count < 2:
                boundary.append(e)
        return boundary

    def get_boundary_vertices(self):
        r"""
            Return
        """
        vertices = []
        for e in self.get_boundary_edges():
            vertices += e
        return list(set(vertices))

    def plot_graph(self):
        r"""
            Plot graph
        """
        plt.clf()
        for e in self.edges:
            v0 = self.vertices[e[0]]
            v1 = self.vertices[e[1]]
            x = [v0[0], v1[0]]
            y = [v0[1], v1[1]]
            plt.plot(x, y, color='k')
        axes = plt.gca()
        vmin = np.amin(self.vertices)
        vmax = np.amax(self.vertices)
        axes.set_xlim([vmin - 0.1, vmax + 0.1])
        axes.set_ylim([vmin - 0.1, vmax + 0.1])
        plt.pause(0.1)