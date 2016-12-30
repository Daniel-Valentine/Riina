"""
Riina - a Python3 node-graph library
-----
Author: Daniel Valentine
E-mail: dval@bu.edu
-----
File  : exceptions.py
Desc  : Custom exception objects
Date  : December 29, 2016
"""

#
# Graph Exceptions
#

class GraphException(Exception):
  pass

class GraphDomainException(GraphException):
  pass

#
# Node Exceptions
#

class NodeException(Exception):
  pass

class NoSuchNodeException(NodeException):
  pass

#
# Edge Exceptions
#

class EdgeException(Exception):
  pass

class NoSuchEdgeException(EdgeException):
  pass

class TooManyEdgesException(EdgeException):
  pass

#
# Path exceptions
#
class InvalidPathException(Exception):
  pass

class NoSuchPathException(InvalidPathException):
  pass

class CycleException(InvalidPathException):
  pass
