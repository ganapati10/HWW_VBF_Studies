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
    
    #Loop under variables
    for variableName, variable in self._variables.iteritems():
      
      tcanvas = ROOT.TCanvas("Signal vs Background_" + variableName, "Signal vs Background_" + variableName, 800, 600) #Iniciate the canvas, one for each variable.
      
      list_canvas [counter] = tcanvas   #list the canvas
      
      counter = counter + 1
      
      #For each cut in cuts.py
      for cutName in self._cuts:
        print "cut_", cutName
        
        if 'cuts' in variable and cutName not in variable['cuts']:
          continue

        if type(fileIn) is not dict and not fileIn.GetDirectory(cutName+"/"+variableName):
          continue
          
        print "variableName =", variableName
            
        # Define variables for store bin content of histograms
        sig  = 0
        bkg  = 0
        data = 0

        
        for samplesName, plotdef in plot:
          
          if 'samples' in variable and sampleName not in variable['samples']:
            continue

          shapeName = cutName+"/"+variableName+'/histo_' + sampleName
          
          if type(fileIn) is dict:
            histo = fileIn[sampleName].Get(shapeName)
          else:
            histo = fileIn.Get(shapeName)
          print ' --> ', histo
          print 'new_histo_' + sampleName + '_' + cutName + '_' + variableName
          histogram = histo.Clone('new_histo_' + sampleName + '_' + cutName + '_' + variableName)
          
          if plotdef['isSignal'] == 1 :
            
            sig = sig + histogram.GetEntries()
            
          elif plotdef['isSignal'] == 0 :
            
            bkg = bkg + histogram.GetEntries()
            
          elif plotdef['isData'] == 1 :
            
            data = data + histogram.GetEntries()
          
        cuts = {}
        
        if os.path.exists(opt.self._cuts) :
          handle = open(opt.self._cuts,'r')
          exec(handle)
          handle.close()      
        
        for cut_k, cut_v in cuts.items():
          if cut_k == cutName:
            txt = cut_v[0].split("&&")
          for cuts_ind in txt:
            cut_ind = cuts_ind.split('>')
            
            if len(var) < 2:
              cut_ind = cuts_ind.split('<')
              
              if len(var) < 2:
                    cut_ind = cuts_ind.split('=')
            
            if cut_ind[0].strip() == cut_k.split("_")[0]:
                
                num = cut_ind[1]
                print(num) ##Cambiar esto, arreglar


         
        
        
            
          
          
      
      
      
      
      
      #Continuara
      
  
  
  
  
  
  
  




















