#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to reconstruct tracking projections from subdirectories and create movies
Requires: readRthRaw.py & theora_png2theora

For each subdir in current directory:
 - for each file matching the given pattern, create a png directory
  - within this directory, reconstruct the projections
  - then create a video of the projections

Note there may be multiple files in the subdirs that match the given pattern.

"""
import os
import glob
import argparse

BASE_DIR=os.path.dirname(__file__)
READRTHRAW_PATH = BASE_DIR + "/readRthRaw.py"
THEORA_EXEC = "theora_png2theora"

def makeMovie(dirname,movie_name):
    dirname = os.path.abspath(dirname)
    callString = THEORA_EXEC + " " + dirname + "/proj%04d.png -f 10 -F 1 -v 10 -o " + movie_name + ".ogv"
    print(callString)
    os.system(callString)

def recon(fname,dirname,ylim):
    startDir = os.getcwd()
    os.chdir(dirname)
    # Run the recon script
    callString = "python " + READRTHRAW_PATH + " " + fname + " -p -s -y " + str(ylim)
    print(callString)
    os.system(callString)
    os.chdir(startDir)

def descend(basedir,filepatt,ylim,movie_prefix):
    if movie_prefix == "":
        movie_prefix = os.path.basename( os.path.abspath(basedir) )
    for lname in os.listdir(basedir):
        subd = os.path.join(basedir, lname)        
        if os.path.isdir(subd):
            print("Subdir: " + subd)
            flist = glob.glob(subd + "/" + filepatt)
            for fname in flist:
                filePrefix = os.path.splitext(fname)[0]
                fileBase = os.path.basename(filePrefix)
                if filePrefix.find("-") != -1:
                    filePrefix = filePrefix[:filePrefix.find("-")]
                    fileBase = fileBase[:fileBase.find("-")]
                pngDir = filePrefix + "-pngDir"
                if not os.path.isdir(pngDir):
                    os.mkdir(pngDir)
                else:
                    print("Warning: will overwrite files in " + pngDir)
                absFname = os.path.abspath(fname)
                recon(absFname,pngDir,ylim)
                subdBase = os.path.basename(subd)
                makeMovie(pngDir,movie_prefix+ "-" + subdBase + "-" + fileBase)

def main():
    parser = argparse.ArgumentParser(description="Reconstruct tracking projections")
    parser.add_argument("--pattern",dest="filepattern",default="cathcoil[4-5]-*.projections",help="file pattern")
    parser.add_argument("--ylim",dest="ylim",type=int,default=800,help="y-axis max")
    parser.add_argument("--prefix",dest="prefix",default="",help="prefix for movie files")
    args = parser.parse_args()
    print("Filepattern: " + args.filepattern)
    descend(".",args.filepattern,args.ylim,args.prefix)
    #makeCats(args.PATTERN1, args.PATTERN2,args.coord)

if __name__ == "__main__":
    main()
