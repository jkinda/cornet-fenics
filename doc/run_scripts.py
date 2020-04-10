#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 15:03:08 2018

@author: Jeremy Bleyer, Ecole des Ponts ParisTech,
Laboratoire Navier (ENPC,IFSTTAR,CNRS UMR 8205)
@email: jeremy.bleyer@enpc.f
"""
import sys
import os

subdir = "demo/"

# Get list of demos (demo name , subdirectory)
demos = [(dI, os.path.join(subdir, dI)) for dI in os.listdir(subdir) if os.path.isdir(subdir)]

# Iterate over demos
for demo, path in demos:
    if os.path.isdir(path):
    # Build list of rst files in demo source directory
        py_files = [f for f in os.listdir(path) if os.path.splitext(f)[1] == ".py" ]
        for f in py_files:
            try:
                ret = os.system("python3 "+os.path.join(path,f))
                if ret == 0:
                    print("Execution of %s went fine." % os.path.join(path,f)) 
                else:
                    raise(ValueError)
            except:
                print("File %s did not run correctly" % os.path.join(path,f))
