#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Read and plot catheter raw data from RthReconImageExporter

    Reads a single raw data file exported from RTHawk. This is assumed to be a catheter
    projection file. The first projection is plotted and subsequent projections can
    be selected for plotting.

    Plots in png format and/or coordinates & physiology data in text format can be output instead using
    the "-p" and "-c" options.

    Note that respiratory data from RTHawk is scaled by 10^5 to give an integer (the respiratory information
    from RTHawk is a float between 0 and 1).

    Usage examples:
        ./readRthRaw.py data/file-0000.projections

        ./readRthRaw.py data/file-0000.projections -p

        ./readRthRaw.py data/file-0000.projections -c

"""
from __future__ import print_function
import sys, struct
import math
import scipy
import scipy.fftpack
from optparse import OptionParser
import os
import snrCalc
import numpy as np
from projPlot import *

if sys.version_info[0] < 3 and sys.version_info[1] < 6:
  raise("Python 2.6+ required...")

float_bytes = 8 #These are being written on a 64-bit system

class RawReader:
    def __init__(self,fname="",float_bytes=8,legacy=False,legacy2=False):
        self.rawFile = fname
        self.legacy = legacy
        self.legacy2 = legacy2
        self.float_bytes = float_bytes
        self.setup()

    def setup(self):
        self.projections = [] # array of tuples, where each tuple is a series of complex floats
        self.projComplex = []
        self.triggerTimes = [] #array of trigger times, one triggerTime per each triplet of projections
        self.respPhases = []
        self.timestamps = []
        self.projNum = 0
        self.projSize = 0
        self.fts = []
        self.projRaw = []     # This will be the same as projComplex but with the x, y and z projections each in their own row
        self.xsize = 0
        self.ysize = 0
        self.zsize = 0
        self.fieldOfView = 0       
 
    def readHeader(self,fp):
      hdr = fp.read(44) # 3 32-bit ints + 1 64-bit int + 3 64-bit floats = 44 bytes
      if not hdr:
        print("reached EOF")
        return (0,0,0,0,0,0,0)
      else:
        xsize = struct.unpack('>i',hdr[0:4])[0]
        ysize = struct.unpack('>i',hdr[4:8])[0]
        zsize = struct.unpack('>i',hdr[8:12])[0]
        fov = struct.unpack('>d',hdr[12:20])[0]
        timestamp = struct.unpack('>q',hdr[20:28])[0]
        trig = struct.unpack('>d',hdr[28:36])[0]
        resp = struct.unpack('>d',hdr[36:44])[0]
        return (xsize,ysize,zsize,fov,trig,resp,timestamp)

    def readLegacy2Header(self,fp):
      hdr = fp.read(36) # 3 32-bit ints + 3 64-bit floats = 36 bytes
      if not hdr:
        print("reached EOF")
        return (0,0,0,0,0,0)
      else:
        xsize = struct.unpack('>i',hdr[0:4])[0]
        ysize = struct.unpack('>i',hdr[4:8])[0]
        zsize = struct.unpack('>i',hdr[8:12])[0]
        fov = struct.unpack('>d',hdr[12:20])[0]
        trig = struct.unpack('>d',hdr[20:28])[0]
        resp = struct.unpack('>d',hdr[28:36])[0]
        return (xsize,ysize,zsize,fov,trig,resp)

    def readLegacyHeader(self,fp):
      hdr = fp.read(20) # 3 32-bit ints + 1 64-bit floats = 20 bytes
      if not hdr:
        print("reached EOF")
        return (0,0,0,0)
      else:
        xsize = struct.unpack('>i',hdr[0:4])[0]
        ysize = struct.unpack('>i',hdr[4:8])[0]
        zsize = struct.unpack('>i',hdr[8:12])[0]
        fov = struct.unpack('>d',hdr[12:20])[0]
        return (xsize,ysize,zsize,fov)
    
    def readFile(self,rawFile = ""):
        if not rawFile:
            rawFile = self.rawFile
        fp = open(rawFile,"rb")
        done = False
        first = True
        self.setup()
        while not done:
          xs = ys = zs = fov = 0
          if self.legacy:
            xs,ys,zs,fov=self.readLegacyHeader(fp)
          elif self.legacy2:
            xs,ys,zs,fov,trig,resp=self.readLegacy2Header(fp)
            self.triggerTimes.append(trig)
            self.respPhases.append(resp)
          else:
            xs,ys,zs,fov,trig,resp,timestamp=self.readHeader(fp)
            self.triggerTimes.append(trig)
            self.respPhases.append(resp)
            self.timestamps.append(timestamp)
          if (xs == 0 or ys == 0 or zs == 0):
            done = True;
            break
          if first:
            self.xsize = xs
            self.ysize = ys
            self.zsize = zs
            self.fieldOfView = fov
          projSize = xs*ys*zs*2
          projByteSize = projSize*float_bytes
          proj = fp.read(projByteSize)
          if proj is None or len(proj) < projByteSize:
            print("Could not read projection " + str(self.projNum) + " stopping here.")
            break
          self.projections.append( struct.unpack('>'+str(projSize)+'d',proj[0:projByteSize]) )
          self.projNum+=1
        print("Read " + str(self.projNum) + " projections...",end='')
        print("x size = " + str(self.xsize) + ", ysize = " + str(self.ysize) 
            + " zsize = " + str(self.zsize) + " fov = " + str(self.fieldOfView))
        # NOTE: each projection in projComplex and projections contains the x, y and z projections
        for proj in range(0,self.projNum):
          self.projComplex.append([])
          for i in range(0,projSize,2):
            self.projComplex[proj].append( complex(self.projections[proj][i],self.projections[proj][i+1]) )

        print("Num projections " + str(len(self.projComplex)))
        for projection in self.projComplex:
          # split into 'ysize' (3) projections
          for y in range(1,self.ysize+1):
            axis = projection[self.xsize*(y-1):self.xsize*y]
            inverseft = scipy.fftpack.ifft(axis) #,npts)
            self.fts.append( scipy.fftpack.fftshift(inverseft) )
            self.projRaw.append(axis)
        print("Num ffts " + str(len(self.fts)))


def getStats(valueList,valueNames):
    """Return stats of given list of sublists, for each index in the sublists

     - Sublists must be of the same length
     - valueNames must be of the same length as sublists
    Example:
    valueList: [ [100,40,20], [110,42,18] ] # [[snr_x, snr_y, snr_z], [snr_x, snr_y, snr_z]]
    valueNames: ["snr_x", "snr_y", "snr_y"]
    """
    means = np.mean(valueList,0)
    mins = np.amin(valueList,0)
    maxs = np.amax(valueList,0)
    stds = np.std(valueList,0)
    statString = ''
    for idx,name in enumerate(valueNames):
        statString += "{0} Mean: {1:.2f} StdDev: {2:.2f} Min: {3:.2f} Max: {4:.2f}\n".format(
        name,means[idx],stds[idx],mins[idx],maxs[idx] )
    return statString
    
def main():
    parser = OptionParser(usage=__doc__)
    parser.add_option("-p", "--plot-save", action="store_true", dest="saveplots",help="save plots to files, no gui", default=False)
    parser.add_option("-c", "--coord-save", action="store_true", dest="savecoords",help="save coordinates to files, no gui", default=False)
    parser.add_option("-l", "--legacy", action="store_true", dest="legacy",help="read legacy files with no trig, resp, or timestamp values", default=False)
    parser.add_option("-m", "--legacy2", action="store_true", dest="legacy2",help="read legacy files with trig and resp but no timestamp values", default=False)
    parser.add_option("-s", "--stemless", action="store_true", dest="stemless", help="stemless - do not display vertical red lines for peak values", default=False)
    parser.add_option("-y", "--ylim", dest="ylim",help="y-axis limit", metavar="YLIM", type="int", default=100)
    parser.add_option("-f", "--statsfile", action="store_true", dest="statsfile",help="save stats to files, no gui", default=False)
    parser.add_option("-x", "--statsprefix",dest="statsPrefix",metavar="PREFIX", help="Prefix for stats file", default="")
    (options,args) = parser.parse_args()

    if (len(args) < 1):
        print(parser.print_help())
        sys.exit(0)
    rawFile = args[0]
    
    rdr = RawReader(rawFile,legacy=options.legacy,legacy2=options.legacy2)
    rdr.readFile()

    if len(rdr.fts) > 0:
      plotter = ProjectionPlot(rdr.fts,rdr.xsize,rdr.fieldOfView,trigTimes=rdr.triggerTimes,
        respArr=rdr.respPhases,timestamps=rdr.timestamps,ysize=rdr.ysize,
        drawStems=(not options.stemless),ylim=options.ylim)

    #Save to files:
    if (options.saveplots or options.savecoords or options.statsfile) and len(rdr.fts) > 0:
      coordFile = None
      allCoords = []
      allSnrs = []
      if options.savecoords:
        fbase,fext = os.path.splitext(rawFile)
        coordFile = open(fbase + '-coords.txt','w')
      for i in range(len(rdr.projComplex)):
        (coords,snrs) = plotter.showProj(i,options.saveplots,options.savecoords,coordFile)
        allCoords.append(coords)
        allSnrs.append(snrs)
        sys.stdout.write("\rSaved projection %i" % i)
        sys.stdout.flush()
      if coordFile:
        coordFile.close()
      print('\n')
      coordStats = getStats(allCoords,["X_coord", "Y_coord", "Z_coord"])
      snrStats = getStats(allSnrs,["SNR_X", "SNR_Y", "SNR_Z"])
      print( coordStats )
      print( snrStats )      
      if options.statsfile:
          fbase,fext = os.path.splitext(rawFile)
          if options.statsPrefix:
              fprefix = options.statsPrefix 
          else:              
              fprefix = os.path.basename(fbase)
          dirname = os.path.dirname(os.path.abspath(rawFile))
          statFile = open(dirname + "/" + fprefix + '-stats.txt','w')
          statFile.write(coordStats + '\n')
          statFile.write(snrStats)
      print("Done.")
    elif len(rdr.fts) > 0: # launch plotter GUI
      plotter.launchGUI()
    else:
      print("Nothing to see here. Exiting.") 

    sys.exit(0)

if __name__ == "__main__":
    main()
