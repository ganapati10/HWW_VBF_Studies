#!/usr/bin/env python

import json
import sys
argv = sys.argv
sys.argv = argv[:1]
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
# Sergio Blanco - 10/12/2020
#


#
# ------------------ Plot of signal vs background ---------------------------------------
#

class Plot_Sig_Bkg:

  _logger = logging.getLogger('Plot_Sig_Bkg')
  
  
  def __init__(self):
    
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
    
    
    
  def makePlot(self, inputFile, outputDirPlots, variables, cuts, samples, plot, legend, groupplot):
    
    print "=================="
    print "==== makePlot ===="
    print "=================="
    
    import LatinoAnalysis.ShapeAnalysis.tdrStyle as tdrStyle
    tdrStyle.setTDRStyle()
        
    ROOT.TGaxis.SetExponentOffset(-0.08, 0.00,"y")
    
    self._variables = variables
    self._samples   = samples
    self._cuts      = cuts

    self._outputDirPlots = outputDirPlots
    os.system ("mkdir " + outputDirPlots + "/") 
    
    ROOT.TH1.SetDefaultSumw2(True)
    
    ROOT.gROOT.cd()
    
    
    fileIn = ROOT.TFile(inputFile, "READ")
    

    #Loop under variables
    for variableName, variable in self._variables.iteritems():
      
      tcanvas = ROOT.TCanvas("Signal vs Background_" + variableName, "Signal vs Background_" + variableName, 800, 600) #Iniciate the canvas, one for each variable.
      
      rang = variable['range']  # The range should be the same as cuts
      tHisto = ROOT.TH1F(variableName, "Signal vs Background at " + variableName, rang[0], rang[1], rang[2])

      
      #For each cut in cuts.py
      for cutName in self._cuts:
        print "cut_", cutName
        
        if 'cuts' in variable and cutName not in variable['cuts']:
          continue

        if type(fileIn) is not dict and not fileIn.GetDirectory(cutName+"/"+variableName):
          continue
        
        #Just compute the variables equal to cuts: variable mjj with cut mjj_1
        if variableName != cutName.strip('_')[0]:      #cutName must be variableName_i
          continue
          
        print "variableName =", variableName
            
        # Define variables for store bin content of histograms
        sig  = 0
        bkg  = 0
        data = 0

        ##-Loop above all samples in plot -> WW, Higgs, DY,..
        for samplesName, plotdef in plot:
          
          if 'samples' in variable and sampleName not in variable['samples']:
            continue

          shapeName = cutName+"/"+variableName+'/histo_' + sampleName
          
          #Check .root file
          if type(fileIn) is dict:
            histo = fileIn[sampleName].Get(shapeName)     #Get the TH1 for each variable, cut and sample.
          else:
            histo = fileIn.Get(shapeName)
          print ' --> ', histo
          print 'new_histo_' + sampleName + '_' + cutName + '_' + variableName
          histogram = histo.Clone('new_histo_' + sampleName + '_' + cutName + '_' + variableName)  #Open the .root file and create histogram
          
          if plotdef['isSignal'] == 1 :
            sig = sig + histogram.GetEntries()
            
          elif plotdef['isSignal'] == 0 :
            bkg = bkg + histogram.GetEntries()
            
          elif plotdef['isData'] == 1 :
            data = data + histogram.GetEntries()
         
        
        ## End of samples loop.     

        tHisto.SetBinContent(cutName.split("_")[1], sig/ROOT.TMath.sqrt(bkg+sig))   #Fill the histograms   
        
        
      #End cuts loop
      
      #TH1F make up
      tHisto.SetMinimum(0.0)
      tHisto.SetMaximum(1.0)
      
      tHisto.SetMarkerStyle(ROOT.kFullCircle)
      tHisto.GetXaxis().SetTitle(variable['xaxis'])
      tHisto.GetYaxis().SetTitle('#frac{S}{\sqrt{B+S}}')
      tHisto.Draw()

      legend = ROOT.TLegend(0.9, 0.87, 0.75, 0.82);
      
      legend.AddEntry(histo, "mjj")
      
      tHisto.SetStats(False)
      
      tHisto.Draw("p")
      
      legend.Draw()
      
      tcanvas.Modified()
      tcanvas.Update()
      tcanvas.Draw() 
      
      tcanvas.SaveAs(self._outputDirPlots + "/" + "Sig_vs_Bkg_" + variableName + ".png")  #Save .png file
    # End variable loop
  #End makePlot function
#End of Plot_Sig_Bkg class  
  
  

#
#
#   PLOT MAKER
#
#   INITIALIZATE THE PLOT CLASS ABOVE
#
#




print '''
--------------------------------------------------------------------------------------------------
   _ \   |         |         \  |         |                
  |   |  |   _ \   __|      |\/ |   _` |  |  /   _ \   __| 
  ___/   |  (   |  |        |   |  (   |    <    __/  |    
 _|     _| \___/  \__|     _|  _| \__,_| _|\_\ \___| _|   
 
--------------------------------------------------------------------------------------------------
'''  

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)

parser.add_option('--inputFile', dest='inputFile', help='input file with histograms', default='input.root')
parser.add_option('--outputDirPlots' , dest='outputDirPlots' , help='output directory', default='./')

hwwtools.addOptions(parser)
hwwtools.loadOptDefaults(parser)
(opt, args) = parser.parse_args()

sys.argv.append( '-b' )
ROOT.gROOT.SetBatch()

print ""
print "          configuration file =", opt.pycfg
print "                        lumi =", opt.lumi
print "                   inputFile =", opt.inputFile
print "              outputDirPlots =", opt.outputDirPlots
print ""


plotter = Plot_Sig_Bkg()

samples = OrderedDict()
if opt.samplesFile == None :
  print " Please provide the samples structure (not strictly needed in mkPlot, since list of samples read from plot.py) "    
elif os.path.exists(opt.samplesFile) :
  handle = open(opt.samplesFile,'r')
  exec(handle)
  handle.close()

cuts = {}
if os.path.exists(opt.cutsFile) :
  handle = open(opt.cutsFile,'r')
  exec(handle)
  handle.close()
  
variables = {}
if os.path.exists(opt.variablesFile) :
  handle = open(opt.variablesFile,'r')
  exec(handle)
  handle.close()

import LatinoAnalysis.ShapeAnalysis.utils as utils

subsamplesmap = utils.flatten_samples(samples)
categoriesmap = utils.flatten_cuts(cuts)

utils.update_variables_with_categories(variables, categoriesmap)

groupPlot = OrderedDict()
plot = {}
legend = {}
if os.path.exists(opt.plotFile) :
  handle = open(opt.plotFile,'r')
  exec(handle)
  handle.close()

plotter.makePlot( opt.inputFile ,opt.outputDirPlots, variables, cuts, samples, plot, legend, groupPlot)

print '.... Now clossing....'
























