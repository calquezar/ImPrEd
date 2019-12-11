#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Force Directed Layout
...description...

 - Data structures:
    * Vertex: List with two elements
    *

REFERENCES:
- :wikipedia:`Force-directed graph drawing`
.. SEEALSO::
    Simonetto, P. , Archambault, D. , Auber, D. and Bourqui, R. (2011),
    ImPrEd: An Improved Force‐Directed Algorithm that Prevents Nodes from Crossing Edges.
    Computer Graphics Forum, 30: 1071-1080. doi:10.1111/j.1467-8659.2011.01956.x
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
import QuadTree as QT
from QuadTree import QPoint, QuadTree
from Graph import SurroundingInfo
from matplotlib import pyplot as plt
from shapely.geometry import LineString, Point, Polygon

class ForceDirectedLayout:

    def __init__(self, graph, beta=1, delta=1, gamma=1, maxIter=20, opt=False):
        r"""
            Constructor
        """
        self.graph = graph
        self.beta = beta
        self.delta = delta
        self.gamma = gamma
        self.maxIter = maxIter
        self.surroundings = self._calculate_surroundings()
        self.opt = opt

    def _calculate_maximum_movements(self, set_of_vertices):
        r"""
            Return the maximum displacements in each direction for each vertex of a given region

            INPUT:

                - vertices

                - region

                - surroundings

            OUTPUT:

                - safe_areas_for_all_vertices
        """
        safe_areas_for_all_vertices = []
        for i in set_of_vertices:
            v = self.graph.vertices[i]
            sv = self.surroundings[i]
            safe_area = []
            for e in sv.edges:
                edge_hindering = self._maximum_displacements(v, e)
                safe_area = self._update_maximum_movements(safe_area, edge_hindering)
            safe_areas_for_all_vertices.append(safe_area)
        return safe_areas_for_all_vertices

    def _calculate_surroundings(self):
        r"""
            Return the surrounding information of each vertex in 'vertices'.

            INPUT:

                - graph
        """
        surroundings = []
        for v in range(len(self.graph.vertices)):
            # Faces
            faces = []
            for r in self.graph.regions:
                if v in r:
                    faces += r
            faces = list(set(faces))  # remove redundancies
            faces.remove(v)  # remove vertex v from the list
            faces.sort()  # sort ascending
            # Connections
            connected_to = []
            for v2 in faces:
                edge = [v, v2] if v < v2 else [v2, v]
                if edge in self.graph.edges:
                    connected_to += [v2]
            connected_to.sort()
            # Surrounding Edges
            boundary = []
            for v1 in range(len(faces) - 1):
                for v2 in range(v1 + 1, len(faces)):
                    if [faces[v1], faces[v2]] in self.graph.edges:  # asumo las aristas se dan con los vertices ordenados de menor a mayor
                        boundary.append([faces[v1], faces[v2]])
            surroundings.append(SurroundingInfo(faces, boundary, connected_to))
        return surroundings

    def _dist(self, v1, v2):
        r"""
            Return the euclidean distance between vertices v1 and v2.
        """
        return math.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2)

    def _forces_calculation(self, set_of_vertices):
        r"""
            Return the corresponding forces for each vertex

            INPUT:

                - vertices

                - region

                - qTree

                - surroundings

                - delta

                - gamma
        """
        centroid = np.mean(self.graph.vertices, axis=0)
        boundary_vertices = self.graph.get_boundary_vertices()
        forces = np.zeros([len(self.graph.vertices), 2])
        for i in range(len(set_of_vertices)):
            v_index = set_of_vertices[i]
            v = np.array(self.graph.vertices[v_index])
            sv = self.surroundings[v_index]
            # fa = np.array([0, 0])
            # fr = np.array([0, 0])
            # fe = np.array([0, 0])
            # fb = np.array([0, 0])
            # Attraction between connected nodes
            for vIdx in sv.connected_to:
                vc = self.graph.vertices[vIdx]
                f = self._node_node_attraction(v, vc)
                # fa[0] += f[0]
                # fa[1] += f[1]
                forces[v_index] += f
            # repulsion between neighbour nodes
            # neighbours = q_tree.findPoints(QPoint(v[0], v[1]), self.delta)
            neighbours = [self.graph.vertices[w] for w in range(len(self.graph.vertices)) if w != i]
            # neighbours = [self.graph.vertices[w] for w in self.surroundings[i].vertices]
            for n in neighbours:
                f = self._node_node_repulsion(v, n)
                # fr[0] += f[0]
                # fr[1] += f[1]
                forces[v_index] += f
            # repulsion with surrounding edges
            for edge in sv.edges:
                edge_coords = [self.graph.vertices[edge[0]], \
                               self.graph.vertices[edge[1]]]
                ve, distance = self.graph._point_edge_projection(v, edge_coords)
                f = self._node_edge_repulsion(v, ve)
                # fe[0] += f[0]
                # fe[1] += f[1]
                forces[edge[0]] -= 0.5*f
                forces[edge[1]] -= 0.5*f
                forces[v_index] += f
            # if v_index in boundary_vertices:
            #     if len(self.surroundings[v_index].connected_to) == 2:
            #         d = self._dist(v, centroid)
            #         scale = self.graph.calculate_scale()
            #         f = (scale/d**2)*(v-centroid)
            #         # print(f)
            #         forces[i] += f

            # attraction to boundary
            # if i in boundary_vertices:
            #     fb = self._boundary_attraction(v)
            # print(fr,fe,fa,fb)
            # calculation of total force
            # total_force = np.array([0, 0])
            # a = 1
            # r = 1
            # e = 1
            # b = 0
            # total_force[0] = a * fa[0] + r * fr[0] + e * fe[0] + b * fb[0]
            # total_force[1] = a * fa[1] + r * fr[1] + e * fe[1] + b * fb[1]
            # print(fa[0], fr[0], fe[0])
            # forces[i] += total_force
        return forces

    def _boundary_attraction(self, v):
        force = [0.0, 0.0]
        center = self.graph.get_center_of_the_graph()
        # calculation of the radius of the circumcicle
        radius = 0
        for i in self.graph.get_boundary_vertices():
            vertex = self.graph.vertices[i]
            dist = self._dist(vertex, center)
            radius = dist if dist > radius else radius
        vectorx = v[0] - center[0]
        vectory = v[1] - center[1]
        module = math.sqrt(vectorx**2+vectory**2)
        force = radius**2-module**2
        angle = math.atan2(vectory, vectorx)
        forcex = force*math.cos(angle)
        forcey = force*math.sin(angle)
        return np.array([forcex, forcey])

    def _node_edge_repulsion(self, v, ve):
        r"""
            Return the repulsion force between vertex v and edge e, given its projection ve.
        """
        d = self._dist(v, ve)
        if d > 0:
            fx = 0
            fy = 0
            # if d < self.gamma:
            fx = (self.gamma/ d**2) * (v[0] - ve[0])
            fy = (self.gamma/ d**2) * (v[1] - ve[1])
                # fx = (((self.gamma - d) ** 2) / d) * (v[0] - ve[0])
                # fy = (((self.gamma - d) ** 2) / d) * (v[1] - ve[1])
            return np.array([fx, fy])
        else:
            epsilon = np.finfo(np.float32).eps
            return np.array([1 / epsilon, 1 / epsilon])

    def _node_node_attraction(self, v1, v2):
        r"""
            Return the attraction force between vertices v1 and v2.
        """
        d = self._dist(v1, v2)
        fx = (d / self.beta) * (v2[0] - v1[0])
        fy = (d / self.beta) * (v2[1] - v1[1])
        return [fx, fy]

    def _node_node_repulsion(self, v1, v2):
        r"""
            Return the repulsion force between vertices v1 and v2.
        """
        d = self._dist(v1, v2)
        if d > 0:
            fx = ((self.delta / d) ** 2) * (v1[0] - v2[0])
            fy = ((self.delta / d) ** 2) * (v1[1] - v2[1])
            return np.array([fx, fy])
        else:
            epsilon = np.finfo(np.float32).eps
            return np.array([1 / epsilon, 1 / epsilon])

    def _maximum_displacements(self, v, e):
        r"""
            Return the maximum displacement (in any direction)
            that vertex 'v' can do without crossing edge 'e'

            INPUT:

                - v:

                - e:

                - vertices: graph

            OUTPUT:

                - safe_displacements:

        """
        edge_coords = [self.graph.vertices[e[0]], self.graph.vertices[e[1]]]  # coordinates of the edge endings
        ve, distance = self.graph._point_edge_projection(v, edge_coords)
        vector_to_edge = self._dist(v, ve)
        vector_direction = (math.atan2(ve[1] - v[1], ve[0] - v[0]) + 2 * math.pi) % (2 * math.pi)  # between [0,2*pi]
        corresponding_sector = int(vector_direction // (math.pi / 4))  # the circumference is divided into 8 sectors
        safe_displacements = []
        for j in range(8):
            diff = (corresponding_sector - j) % 8
            if diff == 0:  # same sector
                sigma = 1
            elif diff == 1 or diff == 2:
                denominator = math.cos(vector_direction - (j + 1) * math.pi / 4)
                if denominator != 0:
                    sigma = 1 / denominator
                else:
                    epsilon = np.finfo(np.float32).eps
                    sigma = 1 / epsilon
            elif diff == 6 or diff == 7:
                denominator = math.cos(vector_direction - j * math.pi / 4)
                if denominator != 0:
                    sigma = 1 / denominator
                else:
                    epsilon = np.finfo(np.float32).eps
                    sigma = 1 / epsilon
            else:
                sigma = 1  # theoretically infinit
            safe_displacements.append(vector_to_edge * sigma / 2)
        return safe_displacements

    def _move(self, set_of_vertices, forces, it):
        r"""
            Move nodes taking into account forces and maximum displacements
        """
        xmin = np.amin(self.graph.vertices[:, 0])
        xmax = np.amax(self.graph.vertices[:, 0])
        ymin = np.amin(self.graph.vertices[:, 1])
        ymax = np.amax(self.graph.vertices[:, 1])
        minimum = min(xmax-xmin, ymax-ymin)
        # self.theta = minimum / (len(self.graph.vertices) ** 2)
        f_square = forces*forces
        total_force = np.sum(f_square, axis=0)
        max_force = np.sqrt(np.max(total_force))
        self.theta = minimum / (max_force*(np.sqrt(len(self.graph.vertices))))
        # print(self.theta)
        colissions = True
        index = 1
        thetas = np.array(len(forces)*[self.theta])
        while colissions:
            # print("Index: " + str(index), self.theta)
            index += 1
            colissions = False
            copy_graph = self.graph.copy()
            for i in range(len(copy_graph.vertices)):
                new_coords = copy_graph.vertices[i] + thetas[i]*forces[i]
                if copy_graph.check_crossings(i, new_coords, self.surroundings):
                    thetas[i] *= 0.5
                    colissions = True
                    # copy_graph.plot_graph(0.1)
                    x = [copy_graph.vertices[i][0], new_coords[0]]
                    y = [copy_graph.vertices[i][1], new_coords[1]]
                    # plt.plot(x, y, color='r')
                    # plt.pause(0.1)
                    # print("colision", i, copy_graph.vertices, new_coords )
                    break
                else:
                    copy_graph.vertices[i] = new_coords
        # print(self.graph.vertices - copy_graph.vertices)
        self.graph = copy_graph


    def _move_nodes(self, set_of_vertices, forces, safe_areas):
        r"""
            Move nodes taking into account forces and maximum displacements
        """
        # self.graph.plot_graph(pause=0.01)
        for i in range(len(set_of_vertices)):
            v = set_of_vertices[i]
            # if v not in self.graph.get_boundary_vertices():
            force_direction = (math.atan2(forces[i][1], forces[i][0]) + 2 * math.pi) % (
                    2 * math.pi)  # angle between [0,2*pi]
            sector = int(force_direction // (math.pi / 4)) % 8
            force_amplitude = math.sqrt(forces[i][0] ** 2 + forces[i][1] ** 2)
            max_displacement = safe_areas[i][sector]
            shift = force_amplitude if abs(force_amplitude) < abs(max_displacement) else max_displacement
            shiftx = shift * math.cos(force_direction)*self.theta
            shifty = shift * math.sin(force_direction)*self.theta
            ########################################################################################
            if not self.opt:
                self.graph.vertices[v][0] += shiftx
                self.graph.vertices[v][1] += shifty
                if self.graph.count_edge_crossings() > 0:  # if the movement adds a crossing ==> revert
                    self.graph.vertices[v][0] -= shiftx
                    self.graph.vertices[v][1] -= shifty
            ########################################################################################
            else:
                new_coords = [0, 0]
                new_coords[0] = self.graph.vertices[v][0] + shiftx
                new_coords[1] = self.graph.vertices[v][1] + shifty
                crossings_found = self.graph.check_crossings(v, new_coords, self.surroundings)
                if not crossings_found: # no crossings
                    self.graph.vertices[v] = new_coords
                    # print("No detecta cruce")
                else:
                    # print("Detecta cruce")
                    None

            ########################################################################################
            # else:
            #     None
                # fx = max_displacement*math.cos(force_direction)*self.theta
                # fy = max_displacement*math.sin(force_direction)*self.theta
                # plt.arrow(self.graph.vertices[v][0] - shiftx,self.graph.vertices[v][1] - shifty, \
                #             fx, fy, width=0.1)
                # plt.pause(0.2)


    # def _point_edge_projection(v, e):
    #     r"""
    #         Return the projection of vertex 'v' into edge 'e'
    #     """
    #     line = LineString([e[0], e[1]])
    #     p = Point(v[0], v[1])
    #     nearest_point = list(line.interpolate(line.project(p)).coords[0])  # using shapely to find the nearest point to line
    #     distance_to_end0 = self._dist(nearest_point, e[0])
    #     distance_to_end1 = self._dist(nearest_point, e[1])
    #     edge_length = self._dist(e[0], e[1])
    #     if distance_to_end0 < edge_length:  # check if the point is inside the segment or not
    #         ve = nearest_point
    #     else:  # get the nearest end
    #         ve = e[0] if (distance_to_end0 < distance_to_end1) else e[1]
    #     return ve

    def run(self, tol=1):
        histGraphs = []
        histSTDs = []
        r"""
            Return the lowest displacements for each direction
        """
        # self.graph.make_convex_boundary()
        for it in range(self.maxIter):
            print("Iter: " + str(it))
            stop = self.graph.get_stop_criteria()
            std, normalized_std, region = self.graph.get_std_area()
            print("std: "+str(normalized_std))
            # print("Iter: " + str(it) + "; Crossings: " + str(self.graph.count_edge_crossings()) + \
            #       "; Regions: " + str(len(self.graph.regions)))
            # q_tree = QuadTree(QT.arrayToList(self.graph.vertices))
            # q_tree.plot()
            # region selection in function of relative area
            # region = self._select_region(self.graph, tol)
            # if not set_of_vertices:  # region == []
            #     print("Finished")
            #     break
            set_of_vertices = range(len(self.graph.vertices))
            # set_of_vertices = region
            # Step 0 Project boundary to circumcircle
            # self.graph.project_boundary_to_circumcircle()
            # Step 1
            forces = self._forces_calculation(set_of_vertices)
            # Step 2
            # safe_displacements = self._calculate_maximum_movements(set_of_vertices)
            # Step 3
            # self._move_nodes(set_of_vertices, forces, safe_displacements)
            self._move(set_of_vertices, forces, it)
            #self.graph.project_boundary_to_circumcircle()
            stop = self.graph.get_stop_criteria()
            # print(stop)
            # crossings = self.graph.count_edge_crossings()
            # if crossings > 0:
            #     break
            # Plot
            if self.graph.plot:
                self.graph.plot_graph(0.01)
                #plt.savefig("Figures/LayoutExamples/100nodes/frame-" + str(it) + ".png", format='png')

            # if self.graph.count_edge_crossings() > 0:
            #     print("Deberia haber avisado antes")
            #     break
            # save current graph
            # histGraphs.append(self.graph.copy())
            histSTDs.append(stop)
        return histSTDs


    def _select_region(self, graph, tol=1):
        r"""
            Return the lowest displacements for each direction
        """
        std, normalized_std, region = graph.get_std_area()
        # print(normalized_std, std)
        if normalized_std < tol:
            region = []
        return region
        # graph.regions.sort(key=graph.abslog_area)  # sort regions by area in ascending order
        # areas = [abs(math.log(graph.area(r))) for r in graph.regions]
        # lowest = areas[0]
        # highest = areas[-1]
        # region = graph.regions[0] if lowest > highest else graph.regions[-1]
        # return region

    def _update_maximum_movements(self, old_limits, new_limits):
        r"""
            Return the lowest displacements for each direction
        """
        if (old_limits == []):
            return new_limits
        else:
            for i in range(len(old_limits)):
                old_limits[i] = old_limits[i] if abs(old_limits[i]) < abs(new_limits[i]) else new_limits[i]
            return old_limits
