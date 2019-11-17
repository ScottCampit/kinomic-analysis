# ODE modeling for transcriptional networks taken from Alon et al.
# Author: Scott Campit
import numpy as np
import pandas as pd
import scipy.integrate
import matplotlib.pyplot as plt

def hill_functions(B, X, K, n):
  active_rate = B*X/(K**n + X**n)
  repress_rate = B/(1+(X/K)**n)

  conc = X/K

  plt.plot(conc, active_rate)
  plot.show()

  plt.plot(conc, repress_rate)
  plt.show()

  return active_rate, repress_rate

B = 40.0
K = 10.0
X = np.linspace(0,2, num=50)
n = 2
hill_functions(B, K, X, n)

########################## Erdos & Renyi Random network ########################
import random

class Node:
    def __init__(self, index):
        self.index = index
        self.neighbors = []

    def __repr__(self):
        return repr(self.index)

def erdos_renyi_rand(n,p):
    """
    The random network is based off Erdos & Renyi's random graph model, where a
    graph is chosen uniformly at random with (n) nodes and (M) edges of a certain
    weight (p).

    Code was adopted from https://github.com/j2kun/erdos-renyi/blob/master/randomgraph.py
    
    """
    vertices = [Node(i) for i in range(n)]
    edges = [(i,j) for i in xrange(n) for j in xrange(i) if random.random() < p]

    for (i,j) in edges:
        vertices[i].neighbors.append(vertices[j])
        vertices[j].neighbors.append(vertices[i])

    return vertices

def dfsComponent(node, visited):
    """
    The depth first search is a recurrent algorithm that picks an unvisited node
    adjacent to the current node. If the node has been visited, we will ignore
    it and visit all adjacent nodes until they all have been visited. Think of
    it as a random walk through a network.
    """
    for v in node.neighbors:
        if v not in visited:
            visted.add(v)
            dfsComponent(v, visited)

def connectedComponents(vertices):
    """
    connectedComponents keeps track of the list of components found so far, and
    the cumulative set of all verticies found with tthe bread first search.
    """
    components = []
    cumulativeVisited = set()

    for v in verticies:
        if v not in cumulativeVisited:
            componentVisited = set([v])
            dfsComponent(v, componentVisited)

            components.append(componentVisited)
            cumulativeVisited |= componentVisited # note that x |= y specifies union

        assert sum(len(c) for c in components) == len(vertices)
        return components

def sizeOfLargestComponent(vertices):
    return max(len(c)) for c in connectedComponents(vertices)

def graphLargestComponentSize(n, theRange):
    """
    graphLargestComponentSize takes the largest component of the independently
    generated graphs as the probability of an edge varies.
    """
    return [(p, sizeOfLargestComponent(randomGraph(n,p))) for p in theRange]

def movingAverage(a, n=3):
    window = np.ones(n) / float(n)
    return np.covolve(a, window, 'same')

def plot(numVertices, xstart, xend, xpts, fil, windowSize=5)
    plt.clf()

    xs = np.linspace(xstart, xend, num=xpts)
    data = graphLargestComponentSize(numVertices, xs)

    ys = [p[1] for p in data]
    movingavg = movingAverage(ys, windowSize)
    newxs = linspace(xstart, xend, num=len(movingavg))

    plt.plot(xs, ys, 'b-', newxs, movingavg, 'r-', lw=2)
    plt.xlabel("p")
    plt.ylabel("Size of largest component")
    plt.ylim(ymax=numVertices)
    plt.xlim((xs[0], xs[-1]))
    plt.title("Phase Transition for Connectivity of a Random Graph")
    plt.savefig(fil)

if __name__ == "__main__":
   plot(50, xstart=0, xend=0.5, xpts=1000, filename="zoomedout-50-1000.png", windowSize=20)
   plot(50, xstart=0, xend=0.15, xpts=1000, filename="zoomedin-50-1000.png", windowSize=20)
   plot(100, xstart=0, xend=0.15, xpts=1000, filename="zoomedin-100-1000.png", windowSize=20)
   plot(500, xstart=0, xend=0.01, xpts=1000, filename="zoomedin-500-1000.png", windowSize=20)

   for i,n in enumerate(range(20, 500, 20), start=1):
      plot(n, xstart=1.0/n, xend=5.0/n, xpts=500, filename="animation/%02d.png" % i, windowSize=20)
