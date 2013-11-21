#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Read and plot catheter raw data from RthReconImageExporter

    Reads a single raw data file exported from RTHawk. This is assumed to be a catheter
    projection file. The first projection is plotted and subsequent projections can
    be selected for plotting.

    Putting another argument (any string) after the first one will enable batch saving of
    plots to png files rather than showing the interactive plotting GUI.
    
    Usage:
        ./readRthRaw.py data/file-0000.projections
    
        readRthRaw.main("data/file-0000.projections")

        ./readRthRaw.py data/file-0000.projections save
        
"""

import sys, struct
import math
import scipy
import scipy.fftpack
import pylab
from matplotlib.widgets import Button

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
        self.index = 0
        self.plots = []
        self.axis = {0:'X',1:'Y',2:'Z'}
        self.stemMarkers = []
        self.stemBase = []
        self.stemLines = []
        
    def showProj(self,frame, saveToFiles=False):
      # fts: all the fourier-transformed projections in one array; x, y, and z each are in their own row
      # frame: which projection to show
      self.index = frame*3
      del self.plots[:]
      self.clearStems()
      #print "Index " + str(index)
      if len(self.fts) < self.index+3 or self.index < 0:
        print "Frame " + str(frame) + " does not exist"
        return
      if saveToFiles:
        pylab.figure(figsize=(13,6))
      for i in range(0,3):
        axes=pylab.subplot('13'+str(1+i))
        pylab.subplots_adjust(bottom=0.2)
        if (self.mode == "phase"):
          self.plots[i].append( pylab.plot(scipy.angle(self.fts[self.index+i])) )
          pylab.title(self.axis[i] + ' Phase Projection');pylab.xticks(self.tick_locs,self.tick_labels)
        else:
          mag = abs(self.fts[self.index+i])
          peak = max(mag)
          peakInd = list(mag).index(peak)
          self.plots.append( pylab.plot(mag) )
          pylab.title(self.axis[i] + ' Magnitude Projection');
          pylab.ylim([0,30]);
          axes.set_autoscaley_on(False);
          pylab.xticks(self.tick_locs,self.tick_labels); 
          stem_marker, stem_lines, stem_base = pylab.stem([peakInd],[peak],'r-','ro');
          self.stemMarkers.append(stem_marker)
          self.stemBase.append(stem_base)
          self.stemLines.append(stem_lines)
          xres = self.fov/self.xsize
          pylab.xlabel(self.axis[i]+':'+'{0:.3}'.format(xres*(peakInd-len(mag)/2))+' mm')
      if saveToFiles:
        pylab.savefig('proj{0:04d}.png'.format(frame))
        self.clearStems()
        pylab.clf()
        pylab.close()
      else:
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

    def clearStems(self):
        for marker in self.stemMarkers:
            marker.remove()
        for base in self.stemBase:
            base.remove()
        for line in self.stemLines:
            line[0].remove()

        del self.stemMarkers[:]
        del self.stemBase[:]
        del self.stemLines[:]

    def redraw(self):
        self.clearStems()
        for i in range(0,3):
            axes = pylab.subplot('13'+str(1+i))
            if (self.mode == "phase"):
                self.plots[i][0].set_ydata(scipy.angle(self.fts[self.index+i]))
            else:
                mag = abs(self.fts[self.index+i])
                self.plots[i][0].set_ydata(mag)
                peak = max(mag)
                peakInd = list(mag).index(peak)
                stem_marker, stem_lines, stem_base = pylab.stem([peakInd],[peak],'r-','ro');
                self.stemMarkers.append(stem_marker)
                self.stemBase.append(stem_base)
                self.stemLines.append(stem_lines)
                xres = self.fov/self.xsize
                pylab.xlabel(self.axis[i]+':'+'{0:.3}'.format(xres*(peakInd-len(mag)/2))+' mm')
            pylab.draw()

    def next(self,event):
        self.index += 3
        self.index = self.index % len(self.fts)
        self.redraw()

    def prev(self,event):
        self.index -= 3
        self.index = self.index % len(self.fts)
        self.redraw()
        
def main(rawFile=None, saveOpt=""):
    if rawFile is None:
        print __doc__
        sys.exit(0)
    fp = open(rawFile,"rb")
    saveToFiles = len(saveOpt) != 0

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
      if proj is None or len(proj) < projByteSize:
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


    if len(fts) > 0:
      plotter = ProjectionPlot(fts,xsize,fieldOfView)

    #Save to files:
    if saveToFiles and len(fts) > 0:
      for i in range(len(projComplex)):
        plotter.showProj(i,True)
        sys.stdout.write("\rSaved projection %i" % i)
        sys.stdout.flush()
      print "\nDone."
    elif len(fts) > 0:
      plotter.showProj(0)
      axprev = pylab.axes([0.7, 0.02, 0.1, 0.075])
      axnext = pylab.axes([0.81, 0.02, 0.1, 0.075])
      bnext = Button(axnext, 'Next')
      bnext.on_clicked(plotter.next)
      bprev = Button(axprev, 'Previous')
      bprev.on_clicked(plotter.prev)

      pylab.show()
    
    sys.exit(0)
    
if __name__ == "__main__":
    main(*sys.argv[1:])