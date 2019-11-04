#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Create plots of tracking sequence projections

- Pass in projections and metadata read from the RTHawk projection file.
- Plots can be saved or displayed in a window.
- Peaks can be displayed or saved into coordinate files; note that peaks are 
calculated using a simple max-peak (not a centroid).
"""

from __future__ import print_function
import pylab
from matplotlib.widgets import Button
import snrCalc
import sys
import math
import numpy as np

if sys.version_info[0] < 3 and sys.version_info[1] < 6:
  raise("Python 2.6+ required...")

class ProjectionPlot:
    def __init__(self,fts,xsize,fov,mode='magnitude',tickDistance=100,
        trigTimes=[],respArr=[], timestamps=[], ysize=3, drawStems=True,
        ylim=100, pgArr=[], ecg1Arr=[], ecg2Arr=[]):
        """ Set up plotter
        fts: all the fourier-transformed projections in one array; x, y, and z each are in their own row
        """
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
        self.trigTimes = trigTimes
        self.useTrig = len(trigTimes) > 0
        self.respArr = respArr
        self.useResp = len(respArr) > 0
        self.timeStamps = timestamps
        self.useTimeStamps = len(timestamps) > 0
        self.ysize = ysize
        self.drawStems = drawStems
        self.ylim = ylim
        self.pgArr = pgArr
        self.ecg1Arr = ecg1Arr
        self.ecg2Arr = ecg2Arr
        self.usePG = len(pgArr) > 0
        self.useECG1 = len(ecg1Arr) > 0
        self.useECG2 = len(ecg2Arr) > 0
        

    def showProj(self,frame, savePlots=False, saveCoords=False, coordFile=None):
      """
      Display / save one projection
      - frame: which projection to show
      """
      self.index = frame*self.ysize
      useTrig = False
      if self.useTrig:
        trig = self.trigTimes[frame]
      if self.useResp:
        resp = self.respArr[frame]
      if self.useTimeStamps:
        timestamp = self.timeStamps[frame]
      del self.plots[:]
      self.clearStems()
      if len(self.fts) < self.index+self.ysize or self.index < 0:
        print("Frame " + str(frame) + " does not exist")
        return
      coords = []
      snrs = []
      if savePlots:
        pylab.figure(figsize=(13,6))
      for i in range(0,self.ysize):
        axes=pylab.subplot('1'+str(self.ysize)+str(1+i))
        pylab.subplots_adjust(bottom=0.2)
        if (self.mode == "phase"):
          self.plots[i].append( pylab.plot(scipy.angle(self.fts[self.index+i])) )
          pylab.title(self.axis[i] + ' Phase Projection');pylab.xticks(self.tick_locs,self.tick_labels)
        else:
          mag = np.abs(self.fts[self.index+i])
          peak = max(mag)
          peakInd = list(mag).index(peak)
          self.plots.append( pylab.plot(mag) )
          snr = snrCalc.getSNR(mag,peak)
          snrs.append(snr)
          pylab.title(self.axis[i] + ' Magnitude Projection');
          pylab.ylim([0,self.ylim]);
          axes.set_autoscaley_on(False);
          pylab.xticks(self.tick_locs,self.tick_labels);
          if self.drawStems:
           stem_marker, stem_lines, stem_base = pylab.stem([peakInd],[peak],'r-','ro');
           pylab.setp(stem_marker,alpha=0.4)
           pylab.setp(stem_lines,alpha=0.4)
           self.stemMarkers.append(stem_marker)
           self.stemBase.append(stem_base)
           self.stemLines.append(stem_lines)
          xres = self.fov/self.xsize
          coords.append(xres*(peakInd-len(mag)/2))
          pylab.xlabel(self.axis[i]+': {0:.2f} mm, SNR: {1}'.format(coords[i],int(round(snr))))
      if savePlots:
        pylab.savefig('proj{0:04d}.png'.format(frame))
        self.clearStems()
        pylab.clf()
        pylab.close()
      if saveCoords and not coordFile is None:
        coordFile.write("%0.1f %0.1f %0.1f %d" % (coords[0], coords[1], coords[2], min(snrs)))
        if self.useTimeStamps:
          coordFile.write(" %d" % (timestamp))
        if self.useTrig:
          coordFile.write(" %d" % (trig))
        if self.useResp:
          coordFile.write(" %d" % (resp * (10**5)))
        if self.usePG:
          coordFile.write(" %f" % (self.pgArr[frame]))
        if self.useECG1:
          coordFile.write(" %f" % (self.ecg1Arr[frame]))
        if self.useECG2:
          coordFile.write(" %f" % (self.ecg2Arr[frame]))
        coordFile.write("\n")
      elif not savePlots and not saveCoords:
        pylab.draw()
      return(coords,snrs)

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
        for i in range(0,self.ysize):
            axes = pylab.subplot('1'+str(self.ysize)+str(1+i))
            if (self.mode == "phase"):
                self.plots[i][0].set_ydata(scipy.angle(self.fts[self.index+i]))
            else:
                mag = np.abs(self.fts[self.index+i])
                self.plots[i][0].set_ydata(mag)
                peak = max(mag)
                peakInd = list(mag).index(peak)
                snr = snrCalc.getSNR(mag,peak)
                if self.drawStems:
                  stem_marker, stem_lines, stem_base = pylab.stem([peakInd],[peak],'r-','ro');
                  pylab.setp(stem_marker,alpha=0.4)
                  pylab.setp(stem_lines,alpha=0.4)
                  self.stemMarkers.append(stem_marker)
                  self.stemBase.append(stem_base)
                  self.stemLines.append(stem_lines)
                xres = self.fov/self.xsize
                pylab.xlabel(self.axis[i]+': {0:.2f} mm, SNR: {1}'.format(
                    xres*(peakInd-len(mag)/2),int(round(snr)) ))
            pylab.draw()

    def next(self,event):
        self.index += 3
        self.index = self.index % len(self.fts)
        self.redraw()

    def prev(self,event):
        self.index -= 3
        self.index = self.index % len(self.fts)
        self.redraw()
    
    def launchGUI(self):
        self.showProj(0)
        axprev = pylab.axes([0.7, 0.02, 0.1, 0.075])
        axnext = pylab.axes([0.81, 0.02, 0.1, 0.075])
        bnext = Button(axnext, 'Next')
        bnext.on_clicked(self.next)
        bprev = Button(axprev, 'Previous')
        bprev.on_clicked(self.prev)
        pylab.show()