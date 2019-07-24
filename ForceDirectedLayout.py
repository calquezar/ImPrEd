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
import QuadTree as QT
from QuadTree import QPoint, QuadTree
from Graph import SurroundingInfo

class ForceDirectedLayout:

    def __init__ (self, graph, delta=1, gamma=1, maxIter=20 ):
        r"""
            Constructor
        """
        self.graph = graph
        self.delta = delta
        self.gamma = gamma
        self.maxIter = maxIter

    def _calculate_maximum_movements(self, vertices, region, surroundings):
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
        for i in region:
            v = vertices[i]
            sv = surroundings[i]
            safe_area = []
            for e in sv.edges:
                edge_hindering = self._maximum_displacements(v, e, vertices)
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

    def _forces_calculation(self, region, qTree, surroundings, delta, gamma):
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
        boundary_vertices = self.graph.get_boundary_vertices()
        forces = []
        for i in region:
            v = self.graph.vertices[i]
            sv = surroundings[i]
            fa = [0, 0]
            fr = [0, 0]
            fe = [0, 0]
            fb = [0, 0]
            # Attraction between connected nodes
            for vIdx in sv.connected_to:
                vc = self.graph.vertices[vIdx]
                f = self._node_node_attraction(v, vc, delta)
                fa[0] += f[0]
                fa[1] += f[1]
            # repulsion between neighbour nodes
            neighbours = qTree.findPoints(QPoint(v[0], v[1]), delta)
            for n in neighbours:
                f = self._node_node_repulsion(v, n.toArray(), delta)
                fr[0] += f[0]
                fr[1] += f[1]
            # repulsion with surrounding edges
            for edge in sv.edges:
                edge_coords = [self.graph.vertices[edge[0]],\
                               self.graph.vertices[edge[1]]]
                ve, distance = self.graph._point_edge_projection(v, edge_coords)
                f = self._node_edge_repulsion(v, ve, gamma)
                fe[0] += f[0]
                fe[1] += f[1]
            # attraction to boundary
            if i in boundary_vertices:
                f = self._boundary_attraction(v)
                fb[0] = f[0]
                fb[1] = f[1]
#            print(fr,fe,fa,fb)
            # calculation of total force
            total_force = [0, 0]
            total_force[0] = 2 * fa[0] + fr[0] + fe[0] + fb[0]
            total_force[1] = 2 * fa[1] + fr[1] + fe[1] + fb[1]
            forces.append(total_force)
        return forces

    def _boundary_attraction(self, v):
        
        boundary_rect = self.graph.get_boundary_rect()
        force = [0.0, 0.0]
        for i in range(4):
            edge = [boundary_rect[i], boundary_rect[(i+1) % 4]]
            ve, dist = self.graph._point_edge_projection(v, edge)
            f = self._node_node_attraction(v, ve, 2)
            force[0] += f[0]
            force[1] += f[1]
        return force
                
        
    def _node_edge_repulsion(self, v, ve, gamma=1):
        r"""
            Return the repulsion force between vertex v and edge e, given its projection ve.
        """
        d = self._dist(v, ve)
        if d > 0:
            fx = 0
            fy = 0
            if d < gamma:
                fx = (((gamma - d) ** 2) / d) * (v[0] - ve[0])
                fy = (((gamma - d) ** 2) / d) * (v[1] - ve[1])
            return [fx, fy]
        else:
            epsilon = np.finfo(np.float32).eps
            return [1 / epsilon, 1 / epsilon]

    def _node_node_attraction(self, v1, v2, delta=1):
        r"""
            Return the attraction force between vertices v1 and v2.
        """
        d = self._dist(v1, v2)
        fx = (d / delta) * (v2[0] - v1[0])
        fy = (d / delta) * (v2[1] - v1[1])
        return [fx, fy]

    def _node_node_repulsion(self, v1, v2, delta=1):
        r"""
            Return the repulsion force between vertices v1 and v2.
        """
        d = self._dist(v1, v2)
        if d > 0:
            fx = ((delta / d) ** 2) * (v1[0] - v2[0])
            fy = ((delta / d) ** 2) * (v1[1] - v2[1])
            return [fx, fy]
        else:
            epsilon = np.finfo(np.float32).eps
            return [1 / epsilon, 1 / epsilon]

    def _maximum_displacements(self, v, e, vertices):
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
        edge_coords = [vertices[e[0]], vertices[e[1]]]  # coordinates of the edge endings
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
                sigma = 20  # theoretically infinit
            safe_displacements.append(vector_to_edge * sigma / 2)
        return safe_displacements

    def _move_nodes(self, graph, region, forces, safe_areas):
        r"""
            Move nodes taking into account forces and maximum displacements
        """
        for i in range(len(region)):
            v = region[i]
            # if v not in graph.get_boundary_vertices():
            force_direction = (math.atan2(forces[i][1], forces[i][0]) + 2 * math.pi) % (2 * math.pi)  # angle between [0,2*pi]
            sector = int(force_direction // (math.pi / 4)) % 8
            force_amplitude = math.sqrt(forces[i][0] ** 2 + forces[i][1] ** 2)
            max_displacement = safe_areas[i][sector]
            shift = force_amplitude if abs(force_amplitude) < abs(max_displacement) else max_displacement
            shiftx = shift * math.cos(force_direction)
            shifty = shift * math.sin(force_direction)
            graph.vertices[v][0] += shiftx
            graph.vertices[v][1] += shifty
            if graph.count_edge_crossings() > 0: #if the movement adds a crossing ==> revert
                graph.vertices[v][0] -= shiftx
                graph.vertices[v][1] -= shifty
            if graph.plot:
                graph.plot_graph()

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
        r"""
            Return the lowest displacements for each direction
        """

        surroundings = self._calculate_surroundings()
        for it in range(self.maxIter):
            qTree = QuadTree(QT.arrayToList(self.graph.vertices))
#            qTree.plot()
            #####################################################
            # region selection in function of relative area
            region = self._select_region(self.graph, tol)
            if region == []:
                print("Finished")
                break
            #####################################################
            # Step 1
            forces = self._forces_calculation(region, qTree, surroundings, self.delta, self.gamma)
            #  #Step 2
            safe_displacements = self._calculate_maximum_movements(self.graph.vertices, \
                                                                                  region, surroundings)
            #  #Step 3
            self._move_nodes(self.graph, region, forces,  safe_displacements)
            # crossings = self.graph.count_edge_crossings()
            # if crossings > 0:
            #     break
            print("Iter: "+str(it)+"; Crossings: " + str(self.graph.count_edge_crossings())+ \
                  "; Regions: " + str(len(self.graph.regions)))


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
