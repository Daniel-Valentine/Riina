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

class NoSuchNodeException(Exception):

  def __init__(self, message = "Node ", errors = {}, *args):
    super().__init__(message, errors, *args)
