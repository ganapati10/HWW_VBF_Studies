#!/usr/bin/env python

import json
import sys
from sys import exit
import ROOT
import optparse
import LatinoAnalysis.Gardener.hwwtools as hwwtools
import os.path
import string
import logging
import LatinoAnalysis.Gardener.odict as odict
import traceback
from array import array
from collections import OrderedDict
import math
import numpy as np
import root_numpy as rnp




#
# ------------------ Plot of signal vs background ---------------------------------------
#

class Plot_Sig_Bkg:

  _logger = logging.getLogger('Plot_Sig_Bkg')
  
  
  def __init__(self):
    
    self._tag = ''
    
    variables = {}
    self._variables = variables
    
    cuts = {}
    self._cuts = cuts
    
    samples = OrderedDict()
    self._samples = samples
    
    
    self._plotLinear = True
    self._plotLog = False

    outputDirPlots = {}
    self._outputDirPlots = outputDirPlots
    
    
    self._fileFormats = ['png']
    
    
  def makePlot(self, inputFile, outputDirPlots, variables, cuts, samples, plot, nuisances, legend, groupplot)
    
    print "=================="
    print "==== makePlot ===="
    print "=================="
    
    self._variables = variables
    self._samples   = samples
    self._cuts      = cuts

    self._outputDirPlots = outputDirPlots
    os.system ("mkdir " + outputDirPlots + "/") 
    
    ROOT.TH1.SetDefaultSumw2(True)
    
    ROOT.gROOT.cd()
    
    
    if os.path.isdir(inputFile):
      # ONLY COMPATIBLE WITH OUTPUTS MERGED TO SAMPLE LEVEL!!
      fileIn = {}
      allFiles = os.listdir(inputFile)
      for sampleName in self._samples:
        fileIn[sampleName] = ROOT.TFile.Open(inputFile+'/plots_%s_ALL_%s.root' % (self._tag, sampleName))
        if not fileIn[sampleName]:
          raise RuntimeError('Input file for sample ' + sampleName + ' missing')
      if os.path.exists(inputFile+'/plots_total.root'):
        fileIn['total'] = ROOT.TFile.Open(inputFile+'/plots_total.root')
          
    else:
      fileIn = ROOT.TFile(inputFile, "READ")
    
    
    
    counter = 0
    
    list_canvas = {}
    
    for variableName, variable in self._variables.iteritems():
      
      tcanvas = ROOT.TCanvas("Signal vs Background_" + variableName, "Signal vs Background_" + variableName, 800, 600)
      
      list_canvas [counter] = tcanvas
      
      
      
      
      
      #Continuara
      
  
  
  
  
  
  
  




















