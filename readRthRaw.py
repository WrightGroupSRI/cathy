#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Read and plot catheter raw data from RthReconImageExporter

    Reads a single raw data file exported from RTHawk. This is assumed to be a catheter
    projection file. The first projection is plotted and subsequent projections can
    be selected for plotting.
    
    Usage:
        ./readRthRaw.py data/file-0000.projections
    
        readRthRaw.main("data/file-0000.projections")
        
"""

import sys, struct
import math
import scipy
import scipy.fftpack
import pylab

float_bytes = 8 #These are being written on a 64-bit system

def readHeader(fp):
  hdr = fp.read(20) # 3 32-bit ints + 1 64-bit float = 20 bytes
  if not hdr:
    print "reached EOF"
    return (0,0,0,0)
  else:
    xsize = struct.unpack('>i',hdr[0:4])[0]
    ysize = struct.unpack('>i',hdr[4:8])[0]
    zsize = struct.unpack('>i',hdr[8:12])[0]
    fov = struct.unpack('>d',hdr[12:20])[0]
    return (xsize,ysize,zsize,fov)

class ProjectionPlot:
    def __init__(self,fts,xsize,fov,mode='magnitude',tickDistance=100):
        self.fts = fts
        self.xsize = xsize
        self.fov = fov
        self.mode = mode
        self.tickDistance = tickDistance
        self.makeTicks()
        
    def showProj(self,frame):
      # fts: all the fourier-transformed projections in one array; x, y, and z each are in their own row
      # frame: which projection to show
      index = frame*3
      #print "Index " + str(index)
      axis = {0:'X',1:'Y',2:'Z'}
      if len(self.fts) < index+3 or index < 0:
        print "Frame " + str(frame) + " does not exist"
        return
    
      for i in range(0,3):
        axes=pylab.subplot('13'+str(1+i))
        if (self.mode == "phase"):
          pylab.plot(scipy.angle(self.fts[index+i]));pylab.title(axis[i] + ' Phase Projection');pylab.xticks(self.tick_locs,self.tick_labels)
        else:
          mag = abs(self.fts[index+i])
          peak = max(mag)
          peakInd = list(mag).index(peak)
          pylab.plot(mag);pylab.title(axis[i] + ' Magnitude Projection');
          pylab.ylim([0,150]);
          axes.set_autoscaley_on(False);
          pylab.xticks(self.tick_locs,self.tick_labels); 
          pylab.stem([peakInd],[peak],'r-','ro');
          xres = self.fov/self.xsize
          pylab.xlabel(axis[i]+':'+'{0:.3}'.format(xres*(peakInd-len(mag)/2))+' mm')
    
      pylab.draw()
    
    def makeTicks(self):
        zeroPixel = (self.xsize+1)/2.0
        tickIncr = self.tickDistance/(self.fov/self.xsize)
        numTicks = int(math.floor(self.xsize / tickIncr))
        self.tick_locs = []
        self.tick_labels = []
        currLabel = -1*self.fov/2.0
        currLoc = 0
        for i in range(0,numTicks):
            self.tick_labels.append(str(currLabel))
            currLabel += self.tickDistance
            self.tick_locs.append(currLoc)
            currLoc += tickIncr
            

def main(rawFile=None):
    if rawFile is None:
        print __doc__
        sys.exit(0)
    fp = open(rawFile,"rb")

    done = False
    projections = [] # array of tuples, where each tuple is a series of complex floats
    projComplex = []
    projNum = 0
    projSize = 0
    first = True
    xsize = ysize = zsize = fieldOfView = 0;
    while not done:
      xs,ys,zs,fov=readHeader(fp)
      if (xs == 0 or ys == 0 or zs == 0):
        done = True;
        break
      if first:
        xsize = xs
        ysize = ys
        zsize = zs
        fieldOfView = fov
      projSize = xs*ys*zs*2
      projByteSize = projSize*float_bytes
      proj = fp.read(projByteSize)
      if not proj:
        print "Could not read projection " + str(projNum) + " stopping here."
        break
      projections.append( struct.unpack('>'+str(projSize)+'d',proj[0:projByteSize]) )
      projNum+=1
    print "Read " + str(projNum) + " projections...",
    print "x size = " + str(xsize) + ", ysize = " + str(ysize) + " zsize = " + str(zsize) + " fov = " + str(fieldOfView)
    # NOTE: each projection in projComplex and projections contains the x, y and z projections
    for proj in xrange(0,projNum):
      projComplex.append([])
      for i in xrange(0,projSize,2):
        projComplex[proj].append( complex(projections[proj][i],projections[proj][i+1]) )

    fts = []
    projRaw = []                # This will be the same as projComplex but with the x, y and z projections each in their own row
    print "Num projections " + str(len(projComplex))
    for projection in projComplex:
      # split into 'ysize' (3) projections
      for y in range(1,ysize+1):
        axis = projection[xsize*(y-1):xsize*y]
        inverseft = scipy.fftpack.ifft(axis) #,npts)
        fts.append( scipy.fftpack.fftshift(inverseft) )
        projRaw.append(axis)
    print "Num ffts " + str(len(fts))

    plotnum = 0
    pylab.ion()
    pylab.figure(plotnum)
    plotter = ProjectionPlot(fts,xsize,fieldOfView)
    plotter.showProj(0)
      
if __name__ == "__main__":
    main(*sys.argv[1:])