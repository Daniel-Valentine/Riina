"""
Riina - a Python3 node-graph library
-----
Author: Daniel Valentine
E-mail: dval@bu.edu
-----
File  : graph.py
Desc  : Basic graph, node, edge implementations
Date  : December 29, 2016
"""

import random as R
import heapq  as H
import math   as M
import sys
import os

from exceptions import *

class Graph(object):

  """
  Graph()                         -> null nodegraph
  Graph(nodes)                    -> nodegraph with collection of unconnected nodes (as ints)
  Graph(nodes, *edges)            -> undirected, possibly cyclic graph with edges given as (int node_1, int node_2[, int weight])
  Graph(nodes, *edges, **kwargs)  -> keyword arguments: visited (directed: bool, source: int, sink: int) 
  """

  # NB only a directed graph!
  
  def __init__(self, V, *E, directed = False, source = -1, sink = -2):
    """ Generate a graph G = (V, E) """

    assert(source != sink)

    self.dir = directed
    self.N = []
    self.S = -1
    self.T = -1

    assert(
      (type(V) == list) &
      (type(E) == list) &
      (len(E) <= len(V)*(len(V)-1)>>1)
    )

    if type(V) != list:
      raise GraphDomainException("The given node set must be a list object")

    # populate node list
    for p in V:
      if type(p) != int:
        raise GraphDomainException("Non-int {} in vertex set".format(p))
      node = self.Node(p)
      if directed:
        if source == p:
          self.S = node
        if sink == p:
          self.T = node
      self.N += [node]
      self.N.sort()

    # populate edge list
    self.g = self.EdgeTable(self, *E)

  def V(self):
    """ Returns vertex set of graph """

    return self.N

  def E(self):
    """ Returns edge set of graph as an EdgeTable object """

    E = []

    for u in g:
      for e in self.g[u]:
        E += [(u, e[0], e[1])]

    return E

  def resetVisits(self):
    """ Unvisit all nodes in the graph """
    for v in self.N:
      v.unvisit()

  def setSource(self, i):
    """ Specify a source node (e.g. for Ford-Fulkerson, Bellman-Ford, or Dijkstra algorithms) """
    s = self.findNode(i)
    assert(type(s) == self.Node)
    self.S = s

  def setSink(self, i):
    """ Specify a sink node (e.g. for Ford-Fulkerson, Bellman-Ford, or Dijkstra algorithms) """
    t = self.findNode(i)
    assert(type(t) == self.Node)
    self.T = t

  def findNode(self, u):
    """ Binary search search for node """

    assert(type(u) == int)
    
    m = 0
    M = len(self.N)-1
    
    while M >= m:
      i = (m+M)>>1
      v = self.N[i].Get()
      if u == v:
        return self.N[i]
      elif u > v:
        m += 1
      else:
        M -= 1

    print("Not found!\n")
    return -1

  def exists(self, u):
    """ Returns True if node u in graph """

    assert(type(u) == int)
    return type(self.findNode(u)) == self.Node

  def reverseEdge(self, u, v):
    """ In a directed graph, reverse the edge u->v """

    assert(type(u) == int)
    assert(type(v) == int)
    assert(self.exists(u))
    assert(self.exists(v))

    u = self.findNode(u)
    v = self.findNode(v)
    weight = -1

    for e in self.g[u]:
      if v == e[0]:
        weight = e[1]
        self.g[u].remove(e)
        if len(self.g[u]) == 0:
          self.g.pop(u)
        break

    self.g.add((v, u, weight))

    return

  class algorithms:

    """
    Graph algorithms
    """

    def FordFulkerson(G):
      return NotImplemented
    
    def findPath(G,                         # the directed graph itself
                 s = -1,                    # origin node; if not specified, will start from defined source node, or raise NoSuchNodeException
                 t = -1,                    # destination node; if not specified, will attempt to reach defined sink node, or raise NoSuchNodeException
                 markPathVisited = False,   # if True, will mark all nodes in discovered path visited
                 paradigm = 0               # edge selection paradigm; default is 0: DFS-based Greedy
                 ):
      """
      Probe a directed graph for a path from node s to node t

      Required arguments:
        * s - origin node; if not specified, will start from defined source node, or raise NoSuchNodeException
        * t - destination node; if not specified, will attempt to reach defined sink node, or raise NoSuchNodeException

      Keyword arguments:
        * markPathVisited - if True, will mark all nodes in discovered path visited
        * paradigm        - edge selection paradigm; available options:
          + 0:  DFS-based Greedy, default
          + 1:  Shortest path (Dijkstra); throws exception if negative edge detected
          + 2:  Shortest path (Bellman-Ford)
          + 3:  DP weight-maximizing, ideal for Ford-Fulkerson
      """

      assert(G.dir), "G is not a directed graph" # for now make sure this only runs if graph is directed

      if s != -1:
        u = G.findNode(s)
        assert(type(u) != int) # raise error if source node does not exist
      else:
        u = G.S
        
      if t != -1:
        t = G.findNode(t)
        assert(type(t) != int) # raise error if destination node does not exist
      else:
        t = G.T
      
      S = [u]
      candidate = [u] # stack of last visited nodes with > 1 edges
      path = []

      while len(S) > 0:
        v = S.pop()
        path += [v]
        #print(path)
        edgeCount = 0

        if v == t and len(path) > 1:
          v.visit()
          if markPathVisited:
            for i in range(len(path)-1):
              G.g.markPathVisited(path[i],path[i+1])
          return path

        # plumb the graph if there are available (unvisited) edges
        ##      if v in self.g and not v.visited:
        ##        v.visit()
        ##        for n in self.g[v]:
        ##          if not n[0].visited or n[0] == u:
        ##            S.append(n[0])
        ##            edgeCount += 1
        ##        if edgeCount > 1:
        ##          candidate.append(v)

        if G.g.contains(v) and not v.visited:
          v.visit()
          for n in G.g[v]:
            if not n[0].visited or n[0] == u:
              S.append(n[0])
              edgeCount += 1
          if edgeCount > 1:
            candidate.append(v)

        # if there are no available outgoing edges, reverse to the last candiate edge
        else:
          v.visit()
          if len(candidate) == 0: #if we've exhaused all candidate nodes, it's clear there's no path to target node
            return None

          # backtrack to last candidate node
          goto = candidate.pop()
          #print("Reversing to", goto)
          goto.unvisit()
          S.append(goto)

          while len(path) > 0 and path.pop() != goto:
            continue

  class Node:
    """
    Node(i[, **kwargs]) -> node with given parameters
    Available keyword arguments:
      * key       - pointer to key that node should hole (defauult value is None) 
      * value     - pointer to data that node should hold (default value is None)
      * visited   - declare node as visited or not (default value is False)
      * nextNode  - pointer next node in the linked list containing nodes with no incoming edges
    """

    def __init__(self,
                 i,                                 # node id
                 key = None,                        # pointer to key the node holds
                 value = None,                      # pointer to data that node holds
                 visited = False,                   # mark node as unvisited by default
                 nextNode = None,                   # for use in creating linked list of nodes with no incoming edges
                 ):
      assert(type(visited) == bool)
      assert(type(i) == int)
      assert(i >= 0)

      self.i = i
      self.visited = visited
      self.next = nextNode

    def __eq__(self, other):
      return self.i == other.Get()

    def __ne__(self, other):
      return self.i != other.Get()

    def __gt__(self, other):
      return self.i > other.Get()

    def __ge__(self, other):
      return self.i >= other.Get()

    def __lt__(self, other):
      return self.i < other.Get()

    def __le__(self, other):
      return self.i <= other.Get()

    def __cmp__(self, other):
      return self.i - other.Get()

    def __hash__(self):
      return hash(self.i)

    def __repr__(self):

      s = '('
      if self.visited:
        return s+'`'+str(self.i)+')'
      else:
        return s+str(self.i)+')'

    def __int__(self):
      return self.i

    def Get(self):
      """ Returns node's key """
      return self.i

    def Set(self, i):
      """ Assign node to new key """
      self.i = i

    def visit(self):
      """ Mark node as visited """
      self.visited = True

    def unvisit(self):
      """ Mark node as unvisited """
      self.visited = False

  class Edge:

    def __init__(self, u, v, weight = None, directed = False, markAsVisited = False):
      assert(type(u) == Graph.Node)
      assert(type(v) == Graph.Node)

      self.u = u
      self.v = v
      if weight is not None:
        self.weighted = True
        self.w = weight
      else:
        self.weighted = False
        self.w = 0
      self.dir = directed
      self.visited = False

    def __repr__(self):
      if self.dir:
        if weighted:
          return self.u + ' --[' + str(self.w) + ']-> ' + self.v
        return self.u + ' -> ' + self.v
      if weighted:
        return self.u + ' --[' + str(self.w) + ']-- ' + self.v
      return self.u + ' -- ' + self.v
      
    def visit(self):
      self.visited = True
    
    def unvisit(self):
      self.visit = False

  class EdgeTable: # edge tuple for u-v: [node v: Node, edgeVisited: bool[, weight: real]]

    def __init__(self,
                 G,                                 # parent graph
                 *edges,                            # list of edges
                 directed = True,                   # is graph directed or not?
                 orderBy = lambda x: x[0],          #
                 hashf = lambda u, v, SIZE: (8831 * min(int(u),int(v))) % SIZE # hashing
                 ):

      self.edgeLengths = -1

      self.G = G
      self.dir = directed
      self.E = {}

      self.add(*edges)

    def __repr__(self):
      #return self.E.__repr__()
      s = ''

      it0 = iter(self.E)
      
      while True:
        u = next(it0)
        s += str(u) + ": "
        it1 = iter(self.E[u])
        while True:
          s += str(next(it1)[0])
          if it1.__length_hint__() == 0:
            break
          s += ", "
        if it0.__length_hint__() == 0:
          break
        s += '\n'

      return s

    def __getitem__(self, u):
      """
      Returns the list of edges to and from u; if graph is directed, returns all outgoing edges from u
      """

      if type(u) != Graph.Node:
        assert(type(u) == int)
        u = self.G.findNode(u)

      if u not in self.E:
        return [-1]
      return self.E[u]

    def contains(self, u):
      return u in self.E

    def pop(self, key):
      """ Remove all edges from node "key" from edge table """

      if key in self.E:
        self.E.pop(key)

    def get(self, u):
      return self.__getitem__(u)

    def add(self, *edges, markAsVisited = False):

      def findNode(u):
        return self.G.findNode(u)

      def quickParse(v, w):

        if w is not None:
          return [[v, markAsVisited, w]]

        return [[v, markAsVisited]]
    
      for e in edges:

        if self.edgeLengths == -1:
          self.edgeLengths = len(e)
        assert(len(e) == self.edgeLengths)
        
        if len(e) == 3:
          assert(e[2] >= 0) # make sure weight is positive
          w = e[2]
        else:
          w = None

        if e[0] <= e[1] or self.dir:
          if (type(e[0]) == type(e[1]) and type(e[0]) == int):
            u = findNode(e[0])
            v = findNode(e[1])
          else:
            u = e[0]
            v = e[1]
        else:
          if (type(e[0]) == type(e[1]) and type(e[0]) == int):
            u = findNode(e[1])
            v = findNode(e[0])
          else:
            u = e[1]
            v = e[0]
        
        if self.dir: # insert all edges in a directed graph
          if u not in self.E:
            self.E[u] = quickParse(v,w)
            continue
          else:
            for f in self.E[u]:
              if f[0] == v: continue
            self.E[u] += quickParse(v,w)
            continue

        else:
          if u not in self.E and v not in self.E:
            self.E[u] = quickParse(v,w)
            continue
          if u in self.E and v not in E:
            for f in self.E[u]:
              if f[0] == v: continue
            self.E[u] += quickParse(v,w)
            continue
          if v in self.E:
            for f in self.E[v]:
              if f[0] == u: continue
            self.E[v] += quickParse(u,w)
            continue

  class Path:
    """
    Iterator over a specific set of Edge objects corresponding to a path
    """

    def __main__(self, u, v, *edges):
      
      if type(u) != Graph.Node or type(v) != Graph.Node:
        raise PathDomainException("Nodes must be Node objects")
      self.u = u
      self.v = v

      self.edges = []
      
      for e in edges:
        if type(e) != Graph.Edge:
          raise PathDomainExcepion("Edges must be Edge objects")
        edges.append(e)

    def __repr__(self):
      return "Path(" + int(u) + ", " + int(v) + ")"
