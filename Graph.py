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
import math
from matplotlib import pyplot as plt
from shapely.geometry import LineString, Point, Polygon
import statistics

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
        # self.project_to_envelope()

    def area(self, region):
        r"""
            Return the area of the region
        """
        points = [self.vertices[v] for v in region]
        pol = Polygon(points)
        return pol.area

    def abslog_area(self, region):
        r"""
            Return the abs(log(area)) of the given region
        """
        return abs(math.log(self.area(region)))

    def count_edge_crossings(self):
        r"""
            Return the number of crossings in the graph
        """
        crossings = 0
        for e1 in range(len(self.edges)-1):
            for e2 in range(e1+1, len(self.edges)):
                edge1 = LineString([self.vertices[self.edges[e1][0]], self.vertices[self.edges[e1][1]]])
                edge2 = LineString([self.vertices[self.edges[e2][0]], self.vertices[self.edges[e2][1]]])
                if edge1.crosses(edge2):
                    crossings += 1
        return crossings

    def get_std_area(self):
        r"""
            Return
        """
        # self.regions.sort(key=self.abslog_area)  # sort regions by area in ascending order

        self.regions.sort(key=self.area)  # sort regions by area in ascending order
        areas = [self.area(r) for r in self.regions]
        self.mean_area = statistics.mean(areas)
        self.std_area = statistics.stdev(areas)
        smallest = abs(areas[0]-self.mean_area)/self.std_area
        biggest = abs(areas[-1]-self.mean_area)/self.std_area
        region = self.regions[0] if smallest > biggest else self.regions[-1]
        self.normalized_std_area = max(smallest, biggest)
        return self.std_area, self.normalized_std_area, region

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

    def get_envelope(self):
        r"""
            Return the smallest square that contains the graph
        """
        vmin = np.amin(self.vertices)
        vmax = np.amax(self.vertices)
        # corners
        bottom_left = [vmin, vmin]
        bottom_right = [vmin, vmax]
        top_left = [vmin, vmax]
        top_right = [vmax, vmax]
        # list of corners
        points = []
        points.append(bottom_left)
        points.append(bottom_right)
        points.append(top_right)
        points.append(top_left)
        return points

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

    def _point_edge_projection(v, e):
        r"""
            Return the projection of vertex 'v' into edge 'e'
        """
        line = LineString([e[0], e[1]])
        p = Point(v)
        nearest_point = line.interpolate(line.project(p))  # using shapely to find the nearest point to line
        distance_to_end0 = nearest_point.distance(Point(e[0]))
        distance_to_end1 = nearest_point.distance(Point(e[1]))
        edge_length = Point(e[0]).distance(Point(e[1]))
        if distance_to_end0 < edge_length:  # check if the point is inside the segment or not
            ve = list(nearest_point.coords[0])
        else:  # get the nearest end
            ve = e[0] if (distance_to_end0 < distance_to_end1) else e[1]
        dist = Point(ve).distance(Point(v))
        return ve, dist

    # def project_to_envelope(self):
    #     r"""
    #         Return the smallest square that contains the graph
    #     """
    #     envelope = self.get_envelope()
    #     boundary_vertices = self.get_boundary_vertices()
    #     for v in boundary_vertices:
    #         coords = self.vertices[v]
    #         min_distance = -1
    #         for i in range(4):
    #             edge_coords = [envelope[i], envelope[(i+1)%4]]
    #             projection, distance = Graph._point_edge_projection(coords, edge_coords)
    #             if min_distance < 0 or min_distance > distance:
    #                 min_distance = distance
    #                 new_coords = projection
    #         self.vertices[v] = new_coords
