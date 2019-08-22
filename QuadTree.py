#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
QuadTree data structure
...description...

REFERENCES:
-

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

from matplotlib import pyplot as plt
from matplotlib import patches
import math
from shapely.geometry import Polygon


def arrayToList(points):
    return [QPoint(p[0], p[1]) for p in points]


def listToArray(points):
    return [[p.x, p.y] for p in points]


def getMinX(points):
    return min([p.x for p in points])


def getMaxX(points):
    return max([p.x for p in points])


def getMinY(points):
    return min([p.y for p in points])


def getMaxY(points):
    return max([p.y for p in points])


###########################################################################
class QPoint:

    def __init__(self, posX, posY):
        self.x = posX
        self.y = posY

    def __str__(self):
        return "(x: " + str(self.x) + ", y: " + str(self.y) + ")"

    def dist(self, p):
        return math.sqrt((self.x - p.x) ** 2 + (self.y - p.y) ** 2)

    def toArray(self):
        return [self.x, self.y]


###########################################################################
class QNode:

    def __init__(self, minx, miny, width, height):
        self.position = QPoint(minx, miny)  # bottom-left corner
        self.width = width
        self.height = height
        self.points = []
        self.children = []
        self.nodeID = "0"
        self.centerOfMass = []

    def __str__(self, tabs=1):
        output = "Node:" + self.nodeID + ", Number of points: " + str(len(self.points)) + \
                 ", Number of children: " + str(len(self.children))
        if (len(self.children) == 0):
            return output
        else:
            children = ""
            for c in self.children:
                children += "\n" + "\t" * tabs + "|+" + c.__str__(tabs + 1)
            return output + children

    def setCenterOfMass(self):
        if (len(self.points) > 0):
            x = ((float)(sum([p.x for p in self.points]))) / len(self.points)
            y = ((float)(sum([p.y for p in self.points]))) / len(self.points)
            self.centerOfMass = QPoint(x, y)

    def getAllCentersOfMasses(self):
        if (len(self.children) == 0):
            return []
        else:
            centers = []
            for c in self.children:
                centers += c.getCenterOfMass()
            centers += [self.centerOfMass]
            return centers

    def getLeaves(self):
        if (len(self.children) == 0):  # I am a leaf
            return [self]
        else:
            leaves = []
            for c in self.children:
                leaves += c.getLeaves()
            return leaves


###########################################################################
class QuadTree:

    def __init__(self, points, threshold=1):
        """
        :rtype: object
        """
        self.threshold = threshold
        self.root = self.generateTree(points)

    def __str__(self):
        if self.root:  # self.root != []
            return self.root.__str__()
        else:
            return "Empty tree"

    def rectangle(self, x, y, width, height):
        """
            (x3,y3) +-------+(x2,y2)
                    |       |
                    |       |
            (x0,y0) +-------+ (x1,y1)
        """
        x0 = x
        y0 = y
        x1 = x + width
        y1 = y
        x2 = x + width
        y2 = y + height
        x3 = x
        y3 = y + height
        return Polygon([(x0, y0), (x1, y1), (x2, y2), (x3, y3)])

    def findPoints(self, p, radius):
        candidates = self.findPointsSquareAtNode(self.root, p, radius)
        neighbours = []
        for q in candidates:
            if (p.dist(q) <= radius):
                neighbours.append(q)
        return neighbours

    def findPointsSquareAtNode(self, node, p, radius):
        points = []
        npos = node.position
        rect1 = self.rectangle(npos.x, npos.y, node.width, node.height)
        rect2 = self.rectangle(p.x - radius, p.y - radius, 2 * radius, 2 * radius)
        if (rect1.intersects(rect2)):
            if (len(node.points) == 0):  # is an empty leaf
                return points
            elif (len(node.children) == 0):  # is a leaf and has points
                neighbours = [q for q in node.points if
                              ((q.dist(p) < radius) & (q.dist(p) > 0))]  # revisar que p sea distinto de q
                points += neighbours
                return points
            else:  # is a node with more than one point
                for c in node.children:
                    points += self.findPointsSquareAtNode(c, p, radius)
                return points
        else:
            return points

    def generateTree(self, points, x=0, y=0, w=0, h=0, nodeID="0"):

        if (len(points) > self.threshold):
            if (nodeID == "0"):  # root node
                minx = getMinX(points)
                miny = getMinY(points)
                width = getMaxX(points) - minx
                height = getMaxY(points) - miny
            else:
                minx = x
                miny = y
                width = w
                height = h
            node = QNode(minx, miny, width, height)
            node.nodeID = nodeID
            node.points = points
            node.setCenterOfMass()
            """
                        +-------+-------+ ^
                        |       |       | |
                        | ch21  | ch22  | |
                        +-------+-------+ (height)
                        |       |       | |
                        | ch11  | ch12  | |
            (minx,miny) +-------+-------+ v
                        <----(width)---->
            """
            #      print(minx,miny,minx+width,miny+height)
            points11 = [p for p in points if (p.x < minx + width / 2) & (p.y < miny + height / 2)]
            ch11 = self.generateTree(points11, minx, miny, width / 2, height / 2, "11")

            points12 = [p for p in points if (p.x >= minx + width / 2) & (p.y < miny + height / 2)]
            ch12 = self.generateTree(points12, minx + width / 2, miny, width / 2, height / 2, "12")

            points21 = [p for p in points if (p.x < minx + width / 2) & (p.y >= miny + height / 2)]
            ch21 = self.generateTree(points21, minx, miny + height / 2, width / 2, height / 2, "21")

            points22 = [p for p in points if (p.x >= minx + width / 2) & (p.y >= miny + height / 2)]
            ch22 = self.generateTree(points22, minx + width / 2, miny + height / 2, width / 2, height / 2, "22")

            node.children = [ch11, ch12, ch21, ch22]
            return node
        else:
            if (nodeID == "0"):  # root node
                if (len(points) > 0):  # there are points inside the root area
                    minx = QPoint.getMinX(points)
                    miny = QPoint.getMinY(points)
                    width = QPoint.getMaxX(points) - minx
                    height = QPoint.getMaxY(points) - miny
                    node = QNode(minx, miny, width, height)
                else:  # root node and no points
                    node = QNode(0, 0, 0, 0)
            else:  # leaf node
                node = QNode(x, y, w, h)
            node.nodeID = nodeID
            node.points = points
            node.setCenterOfMass()
            return node

    def plot(self):
        fig = plt.figure(figsize=(12, 12))
        plt.title("Quadtree")
        ax = fig.add_subplot(111)
        leaves = self.root.getLeaves()
        for region in leaves:
            ax.add_patch(patches.Rectangle((region.position.x, region.position.y), \
                                           region.width, region.height, fill=False))
        x = [point.x for point in self.root.points]
        y = [point.y for point in self.root.points]
        plt.plot(x, y, 'ro', color='r', markersize=3)
        plt.show()
        return ax
##############################################################
# import random
# random.seed(10)
# points = [QPoint(random.uniform(0, 10), random.uniform(0, 10)) for x in range(10)]
# q = QuadTree(points,1)
##print(q)
#
## test points in a circle
# ax = q.plot()
# x = 8.6
# y = 6
# radius = 4
# print(len(q.findPoints(QPoint(x,y),radius)))
# circle = plt.Circle((x,y), radius, color='r',fill=False)
# ax.add_artist(circle)