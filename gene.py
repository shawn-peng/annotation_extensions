# -*- coding: utf-8 -*-
"""
Created on Sat Jan  5 13:48:49 2019

@author: Shawn
"""

import re
from collections import defaultdict


__all__ = ['Gene']

class Gene:
    def __init__(self):
        self.annotations = {}
        
    def add_annotation(self, gaf_rec):
        ann = Annotation(gaf_rec["GO_ID"])
        ann.add_evidence(gaf_rec["Evidence"])
        for extstr in re.split(",|\|", gaf_rec["Annotation_Extension"]):
            ext = Extension(extstr)
            ann.add_extension(ext)
        
        self.annotations[ann.term] = ann
    
    
    

class Annotation:
    def __init__(self, term):
        self.term = term
        self.evidences = set()
        self.extensions = defaultdict(set)
    
    def add_evidence(self, e):
        self.evidences.add(e)
    
    def add_extension(self, ext):
        self.extensions[ext.rel].add(ext.target)

ext_re = re.compile("(\w+)\((.+)\)")

class Extension(object):
    def __init__(self, s):
        self.rel = ext_re.match(s).group(1)
        self.target = ext_re.match(s).group(2)