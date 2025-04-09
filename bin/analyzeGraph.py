import sys

from matplotlib import pyplot as plt
import numpy as np

import igraph as ig

def local_clustering_coefficient_directed(graph, vertex):
    neighbors_in = set(graph.neighbors(vertex, mode="in"))
    neighbors_out = set(graph.neighbors(vertex, mode="out"))
    neighbors = neighbors_in.union(neighbors_out)
    
    if len(neighbors) < 2:
        return 0.0
    
    possible_edges = len(neighbors) * (len(neighbors) - 1)
    actual_edges = 0
    
    for neighbor in neighbors:
        for other_neighbor in neighbors:
            if neighbor != other_neighbor and (graph.are_adjacent(neighbor, other_neighbor)):
                actual_edges += 1
    
    return actual_edges / possible_edges

gmlfile = sys.argv[1]

graph = ig.Graph.Read_GML(gmlfile)

ig.summary(graph)

print("The graph has " + str(len(graph.vs)) + " vertices")
print("The graph has " + str(len(graph.es)) + " edges")

num_components = len(ig.Graph.connected_components(graph, mode = "weak"))
print("The graph has " + str(num_components) + " weak connected components")

avg_path_length = ig.Graph.average_path_length(graph)
print("Average path length: " + str(avg_path_length))

radius = ig.Graph.radius(graph)
diameter = ig.Graph.diameter(graph)

print("graph radius: " + str(radius))
print("graph diameter: " + str(diameter))

local_clustering_coefficients = [local_clustering_coefficient_directed(graph, v) for v in range(graph.vcount())]
average_local_clustering_coefficient = sum(local_clustering_coefficients) / len(local_clustering_coefficients)
print("graph clustering coefficient: " + str(average_local_clustering_coefficient))

layout = graph.layout("drl")

fig, ax = plt.subplots()
visual_style = {}
visual_style["vertex_size"] = 10
visual_style["label_size"] = 2
visual_style["arrow_width"] = 1

ig.plot(graph, layout=layout, target=ax, **visual_style)

plt.show()
indegrees = []
outdegrees = []

for v in graph.vs:
    indegrees.append(1 + v.indegree())
    outdegrees.append(1 + v.outdegree())

fig, axs = plt.subplots(2, 1, figsize=(16, 10))
fig.suptitle("degree distribution", fontsize=16)
plt.subplot(2, 1, 1)
plt.hist(indegrees, bins=np.logspace(start=np.log10(1), stop=np.log10(1000), num=20))
axs[0].set_ylabel("# observations")
axs[0].set_xlabel("degree")
axs[0].set_xscale("log")
axs[0].set_yscale("log")
axs[0].set_title("indegree")

plt.subplot(2, 1, 2)
plt.hist(outdegrees, bins=np.logspace(start=np.log10(1), stop=np.log10(1000), num=20))
axs[1].set_ylabel("# observations")
axs[1].set_xlabel("degree")
axs[1].set_xscale("log")
axs[1].set_yscale("log")
axs[1].set_title("outdegree")
#hist = ig.Graph.degree_distribution(graph, bin_width=1)


plt.show()
