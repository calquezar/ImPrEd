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
import copy
from matplotlib import pyplot as plt
from shapely.geometry import LineString, Point, Polygon
from shapely.geometry.polygon import LinearRing
from descartes.patch import PolygonPatch
import statistics

class SurroundingInfo:
    def __init__(self, v, e, c):
        self.vertices = v.copy()  # vertices in the boundary of the surrounding area
        self.edges = e.copy()  # edges forming the boundary of the surrounding area
        self.connected_to = c.copy()  # edges that connect the vertex to others in the boundary
#########################################################

def addLimitPoints(points):
    r"""
        Add points in order to close voronoi diagram
    """
    means = np.mean(points, axis=0)
    mins = np.min(points,axis=0)
    maxs = np.max(points,axis=0)
    beta = 2
#    points = np.append(points, [means[0], beta*np.abs(maxs[1])]) # top
#    points = np.append(points, [means[0], -beta*abs(mins[1])]) # bottom
#    points = np.append(points, [-beta*abs(mins[0]), means[1]]) # left
#    points = np.append(points, [beta*abs(maxs[0]), means[1]]) # right
    
    points = np.append(points, [-beta*abs(mins[0]+1), beta*np.abs(maxs[1]+1)]) # top-left
    points = np.append(points, [-beta*abs(mins[0]+1), -beta*abs(mins[1]+1)]) # bottom-left
    points = np.append(points, [beta*abs(maxs[0]+1),beta*np.abs(maxs[1]+1)]) # top-right
    points = np.append(points, [beta*abs(maxs[0]+1), -beta*abs(mins[1]+1)]) # bottom-right
    L = int(points.size)
    points = points.reshape(int(L/2),2)
    return points
    
class Graph:

    def __init__(self, vor, clockwise=False):
        r"""
            Constructor
        """
        self.vertices = vor.vertices
        self.edges = [x for x in vor.ridge_vertices if -1 not in x]
        self.regions = [x for x in vor.regions if -1 not in x]
        self.regions = [x for x in self.regions if x]  # remove empty list
        self.sort_all_regions(clockwise)
        self.plot = False
        # self.project_to_envelope()
        
    def copy(self):
        return copy.deepcopy(self)

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
                    # print(self.edges[e1],self.edges[e2])
        return crossings

    def _dist(self, v1, v2):
        r"""
            Return the euclidean distance between vertices v1 and v2.
        """
        return math.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2)

    def get_boundary_rect(self):
        
        mins = np.min(self.vertices, axis=0)
        maxs = np.max(self.vertices, axis=0)
        rect = []
        rect += [[mins[0], mins[1]]] # bottom-left
        rect += [[maxs[0], mins[1]]] # bottom-right
        rect += [[maxs[0], maxs[1]]] # top-right
        rect += [[mins[0], maxs[1]]] # top-left
        return rect

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
                if self.is_in_region(e, r):
                    count += 1
                if count > 1:
                    break
            if count < 2:
                boundary.append(e)
        return boundary

    def get_boundary_vertices(self, clockwise=False):
        r"""
            Return
        """
        vertices = []
        for e in self.get_boundary_edges():
            vertices += e
        boundary = list(set(vertices)) # remove redundancies
        boundary = self.sort_region_vertices(boundary, clockwise)
        return boundary

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

    def is_in_region(self, edge, region):
        r"""
            Return true if the edge is in the boundary of the region
        """
        isIn = True
        for v in edge:
            if v not in region:
                isIn = False
                break
        return isIn

    def get_region_boundary(self, region):
        edges = []
        for i in range(len(region)-1):
            vertex1 = region[i]
            for j in range(i+1, len(region)):
                vertex2 = region[j]
                edge = [vertex1, vertex2] if vertex1 < vertex2 else [vertex2, vertex1]
                if edge in self.edges:
                    edges.append(edge)
        return edges

    # def get_region_boundary(self, region):
    #     edges = []
    #     for i in range(len(region)):
    #         vertex1 = region[i]
    #         vertex2 = region[(i+1)%len(region)]
    #         edge = [vertex1, vertex2] if vertex1 < vertex2 else [vertex2, vertex1]
    #         if edge in self.edges:
    #             edges.append(edge)
    #     return edges

    def get_region_coords(self, region):

        coords = []
        for v in region:
            coords.append(self.vertices[v])
        return coords

    def is_ccw(self, region):
        coords = self.get_region_coords(region)
        ring = LinearRing(coords)
        if ring.is_ccw:
            return True
        else:
            return False

    def get_regions_by_edge(self, e):

        if e[1] < e[0]:
            e.reverse()
        regions = []
        for r in self.regions:
            boundary = self.get_region_boundary(r)
            if e in boundary:
                regions.append(r)
                if len(regions) == 2:
                    break
        return regions

    def get_consecutive_edge(self, edge, common_vertex, list_edges):

        for e in list_edges:
            if e != edge and common_vertex in e:
                return e

    # def get_edges_by_vertex_in_region(self, vertex, region):
    #
    #     boundary = self.get_region_boundary(region)
    #     edges = []
    #     for e in boundary:
    #         if vertex in e:
    #             edges.append(e)
    #     return edges

    def remove_region(self, region):

        new_graph = self.copy()
        if region in new_graph.regions:
            graph_boundary = new_graph.get_boundary_edges()
            region_boundary = new_graph.get_region_boundary(region)
            for e in region_boundary:
                if e in graph_boundary:
                    new_graph.edges.remove(e)
            new_graph.regions.remove(region)
        return new_graph

    def sort_region_vertices(self, region, clockwise=False):

        # firts get the region boundary
        edges = self.get_region_boundary(region)
        # select the first edge of the list
        e = edges[0]
        # create the new list (initially empty)
        new_ordering = []
        # the first vertex will be the first end of the edge 'e'
        new_ordering.append(e[0])
        # now follow the boundary of the region
        next_vertex = e[1]
        next_edge = self.get_consecutive_edge(e, next_vertex, edges)
        while next_vertex != new_ordering[0]:
            # add new vertex to list
            new_ordering.append(next_vertex)
            # get the other ending of the edge
            next_vertex = next_edge[0] if next_edge[0] != next_vertex else next_edge[1]
            next_edge = self.get_consecutive_edge(next_edge, next_vertex, edges)

        if clockwise:
            if self.is_ccw(new_ordering):
                new_ordering.reverse()
        else:
            if not self.is_ccw(new_ordering):
                new_ordering.reverse()

        return new_ordering

    def sort_all_regions(self, clockwise=False):
        for r in range(len(self.regions)):
            region = self.regions[r]
            self.regions[r] = self.sort_region_vertices(region, clockwise)

    def colour_region(self, region):

        points = []
        for v in region:
            points.append(self.vertices[v])

        polygon = Polygon(points)
        patch = PolygonPatch(polygon, facecolor=[0, 0, 0.5], edgecolor=[0, 0, 0], alpha=0.7, zorder=2)
        ax = plt.gca()
        ax.add_patch(patch)
        plt.pause(0.5)

    def colour_region_edges(self, region, sorted=False, clockwise=False):
        if sorted:
            region = self.sort_region_vertices(region, clockwise)
            for i in range(len(region)):
                e = [region[i], region[(i + 1) % (len(region))]]
                self.colour_edge(e)
        else:
            edges = self.get_region_boundary(region)
            for e in edges:
                self.colour_edge(e)

    def colour_edge(self, e):
        v0 = self.vertices[e[0]]
        v1 = self.vertices[e[1]]
        x = [v0[0], v1[0]]
        y = [v0[1], v1[1]]
        plt.plot(x, y, color='r')
        plt.pause(1)

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

##############################################################################################

    def _point_edge_projection(self, v, e):
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

    def project_boundary_to_circumcircle(self):
        boundary_vertices = self.get_boundary_vertices(clockwise=True)
        N = len(boundary_vertices)  # number or boundary vertices
        # calculation of the center of the circumcicle
        mass_center = [0, 0]
        for v in boundary_vertices:
            mass_center[0] += self.vertices[v][0]/N
            mass_center[1] += self.vertices[v][1]/N
        # center the graph to the origin
        for i in range(len(self.vertices)):
            self.vertices[i][0] -= mass_center[0]
            self.vertices[i][1] -= mass_center[1]
        mass_center = [0, 0]  # Now the center of mass is the origin
        # calculation of the radius of the circumcicle
        radius = 0
        for i in boundary_vertices:
            v = self.vertices[i]
            dist = self._dist(v, mass_center)
            radius = dist if dist > radius else radius
        # calculate the new boundary coordinates
        for i in boundary_vertices:
            v = self.vertices[i]
            angle = math.atan2(v[1], v[0])
            self.vertices[i] = [radius*math.cos(angle), radius*math.sin(angle)]
        # normalize all vertices to the circumference of radius 1
        for i in range(len(self.vertices)):
            self.vertices[i][0] /= radius
            self.vertices[i][1] /= radius

        # for v in boundary_vertices:
        #     self.vertices[v] = new_vertices.pop()

        # self.plot_graph()
        # c = plt.Circle(mass_center, radius)
        # ax = plt.gca()
        # ax.add_artist(c)
