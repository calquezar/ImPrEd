import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from Graph import Graph, addLimitPoints
from ForceDirectedLayout import ForceDirectedLayout
import time
import copy
##########################################################
# # caso importante
# np.random.seed(10)
# case = 1
# if case == 0:
#     points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
#                        [2, 1], [2, 2]])
# elif case == 1: # with 20,2 we found problems
#     points = np.random.random((20, 2))
#
#
# points = addLimitPoints(points)
# vor = Voronoi(points)
# g = Graph(vor)
#############################################################
# caso importante
np.random.seed(10)
case = 1
if case == 0:
    points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
                       [2, 1], [2, 2]])
elif case == 1: # with 20,2 we found problems
    points = np.random.random((20, 2))
elif case == 2:
    points = np.random.random((20, 2))


points = addLimitPoints(points)
vor = Voronoi(points)
g = Graph(vor)

# g.plot_graph()
g.plot = False
# Force algorithm
maxIter = 20
tol = 0.2
beta = g.calculate_scale()  # node_node_attraction
delta = 0.01  # node_node_repulsion
gamma = 0.1  # node_edge_repulsion
f = ForceDirectedLayout(g, beta, delta, gamma, maxIter, opt=True)
f.run()
g = f.graph
#############################################################
# # check region-relative colouring
# for r in g.regions:
#     g.plot_graph()
#     g.colour_region(r)
#     g.colour_region_edges(r, sorted=True, clockwise=False)

# check edge-relative colouring
# for edge in g.edges:
#     g.plot_graph()
#     g.colour_edge(edge)
#     regions = g.get_regions_by_edge(edge)
#     for r in regions:
#         g.colour_region(r)

# boundary checking
# boundary_vertices = g.get_boundary_vertices(clockwise=False)
# for i in range(len(boundary_vertices)):
#     edge = [boundary_vertices[i], boundary_vertices[(i+1)%len(boundary_vertices)]]
#     g.colour_edge(edge)

# # check that regions have been sorted in the graph constructor
# # plotting consecutive vertices in the region list
# for r in range(len(g.regions)):
#     region = g.regions[r]
#     for i in range(len(region)):
#         edge = [region[i], region[(i+1)%len(region)]]
#         g.colour_edge(edge)

# # plot boundary regions only
# boundary_vertices = g.get_boundary_vertices(clockwise=False)
# g.plot_graph()
# for i in range(len(boundary_vertices)):
#     v0 = boundary_vertices[i]
#     v1 = boundary_vertices[(i+1)%len(boundary_vertices)]
#     edge0 = [v0, v1] if v0 < v1 else [v1, v0]
#     region = g.get_regions_by_edge(edge0)[0]
#     g.colour_region(region)


# # # remove boundary regions only, one-by-one
# boundary_vertices = g.get_boundary_vertices(clockwise=False)
# g2 = g.copy()
# for i in range(len(boundary_vertices)):
#     v0 = boundary_vertices[i]
#     v1 = boundary_vertices[(i+1)%len(boundary_vertices)]
#     e = [v0, v1] if v0 < v1 else [v1, v0]
#     if e in g2.get_boundary_edges():
#         region = g2.get_regions_by_edge(e)[0]
#         g2.plot_graph()
#         g2.colour_region(region)
#         g2 = g2.remove_region(region)
#         g2.plot_graph()

g.sort_all_regions(clockwise=False)
def find_basis(g, edge, vertex_ending, generators=[], source_path=[], clockwise=False, initialGraph = None,count=0):

    regions = g.get_regions_by_edge(edge)
    if regions:
        region = regions[0]  # only one region because the edge is on the boundary of the graph
        plt.cla()
        initialGraph.plot_graph(pause=0.1, color=(0, 0, 0, 0.2))
        g.plot_graph()
        # plt.savefig("Figures/frame-" + str(count) + ".png", format='png')
        count += 1
        # g.colour_edge(edge)
        ################################################################
        startv = edge[0] if edge[0] != vertex_ending else edge[1]
        index = region.index(startv)
        # index0 = region.index(edge[0])
        # index1 = region.index(edge[1])
        # if index0 == 0 and index1 == len(region)-1:
        #     shift = index1
        # elif index0 == len(region)-1 and index1 == 0:
        #     shift = index0
        # else:
        #     shift = index0 if index0 < index1 else index1
        shift = index
        region = np.roll(region, -shift).tolist()
        ################################################################
        # closed path around the region of interest
        closedPathRegion = []
        for w in range(len(region)):
            w1 = region[w]
            w2 = region[(w+1)%len(region)]
            e = [w1, w2]
            closedPathRegion += [e]
            # closedPathRegion += [[w1, w2]]
        # geometric generator of the corresponding region from base point
        generator = copy.deepcopy(source_path)
        generator += closedPathRegion
        reversedPath = copy.deepcopy(source_path)
        reversedPath.reverse()
        for e in reversedPath:
            e.reverse()
        generator += reversedPath

        for e in source_path:
            g.colour_edge(e, pause=0.1)
        # plt.savefig("Figures/frame-" + str(count) + ".png", format='png')
        count += 1

        g.colour_region(region)

        for e in closedPathRegion:
            g.colour_edge(e, pause=0.1)
            # plt.savefig("Figures/frame-" + str(count) + ".png", format='png')
            count += 1

        # for e in reversedPath:
        #     g.colour_edge(e, pause=0.1)

        # plt.pause(5)
        # add generator to the list of generators
        generators += [generator]
        # region = region[shift:]+region[:shift]

        # g.plot_graph()
        # for e in source_path:
        #     g.colour_edge(e)
        # for i in range(len(region)):
        #     e = [region[i], region[(i+1)%len(region)]]
        #     e = [e[0], e[1]] if e[0] < e[1] else [e[1], e[0]]
        #     g.colour_edge(e)
        source_path += [[startv, vertex_ending]]
        boundary_edges = g.get_boundary_edges()
        next_edge = g.get_consecutive_edge(edge, vertex_ending, boundary_edges)
        vertex_ending = next_edge[0] if next_edge[0] != vertex_ending else next_edge[1]
        while (g.is_in_region(next_edge, region) or len(g.get_regions_by_edge(next_edge)) == 0) and len(g.regions) > 1:
            source_path += [next_edge]
            if len(g.get_regions_by_edge(next_edge)) == 0:  # edge does not belong to any region
                g.remove_edge(next_edge)
            next_edge = g.get_consecutive_edge(next_edge, vertex_ending, boundary_edges)
            vertex_ending = next_edge[0] if next_edge[0] != vertex_ending else next_edge[1]

        region = np.roll(region, shift).tolist()  # restore the default value of the region
        find_basis(g.remove_region(region, source_path[-1]), next_edge, vertex_ending, generators=generators, source_path=source_path, clockwise=clockwise, initialGraph=initialGraph, count=count)
    else:
        print("No regions")
    return generators

def combineGenerators (g1, g2):
    g1copy = copy.deepcopy(g1)
    g2copy = copy.deepcopy(g2)
    limit = min(len(g1copy), len(g2copy))
    g1copy.reverse()
    i = 0
    while i < limit:
        e1 = g1copy[i]
        e2 = g2copy[i]
        # e2.reverse()
        if e1[0] in e2 and e1[1] in e2:
            # e2.reverse()
            g1copy.remove(e1)
            g2copy.remove(e2)
            i = 0
            limit = min(len(g1copy), len(g2copy))
        else:
            print(e1, e2)
            # e2.reverse()
            break
    g1copy.reverse()
    combined = g1copy + g2copy
    return combined



boundary_vertices = g.get_boundary_vertices(clockwise=False)
v0 = boundary_vertices[0]
v1 = boundary_vertices[1]
edge = [v0, v1] if v0 < v1 else [v1, v0]
basis = find_basis(g, edge, v1, initialGraph=g)
############################################################3
# count = 0
# for generator in basis:
#     g.plot_graph(pause=0.1, color=(0, 0, 0, 0.3))
#     #plt.savefig("Figures/Basis/20nodes/frame-"+str(count)+".png", format='png')
#     count += 1
#     # print(generator)
#     for edge in generator:
#         e = copy.deepcopy(edge)
#         if e[0] > e[1]:
#             e.reverse()
#         g.colour_edge(e, pause=0.1)
#         #plt.savefig("Figures/Basis/20nodes/frame-" + str(count) + ".png", format='png')
#         count += 1
############################################################3
# comb = basis[0]
# count = 0
# for i in range(1, len(basis)):
#     comb = combineGenerators(basis[i], comb)
#     g.plot_graph(pause=0.1, color=(0, 0, 0, 0.3))
#     plt.savefig("Figures/Basis/20NodesBasisComposition/frame-"+str(count)+".png", format='png')
#     count += 1
#     for edge in comb:
#         e = copy.deepcopy(edge)
#         if e[0] > e[1]:
#             e.reverse()
#         g.colour_edge(e, pause=0.1)
#         plt.savefig("Figures/Basis/20NodesBasisComposition/frame-" + str(count) + ".png", format='png')
#         count += 1