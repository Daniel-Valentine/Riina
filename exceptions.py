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
  """ Generic exception for Graph objects """
  pass

class GraphDomainException(GraphException):
  """ Invalid arguments passed into a Graph constructor """
  pass

#
# Node Exceptions
#

class NodeException(Exception):
  """ Generic exception for Node objects """
  pass

class NoSuchNodeException(NodeException):
  """ Node does not exist in graph """
  pass

#
# Edge Exceptions
#

class EdgeException(Exception):
  """ Generic exception for Edge and EdgeTable objects """
  pass

class EdgeDomainException(EdgeException):
  pass

class NoSuchEdgeException(EdgeException):
  pass

class TooManyEdgesException(EdgeException):
  pass

#
# Path exceptions
#

class PathException(Exception):
  """ Generic exception for Path objects """
  pass

class InvalidPathException(PathException):
  """ Path contains a structural problem """
  pass

class PathDomainException(PathException):
  """ Invalid arguments passed into a Path constructor """
  pass

class NoSuchPathException(PathException):
  """ Path does not exist in Graph object """
  pass

class CycleException(PathException):
  """ Path is a cycle """
  pass

class PathBranchException(InvalidPathException):
  """ Attempt to add edge to Path object caused a branch to form in the path """
  pass
