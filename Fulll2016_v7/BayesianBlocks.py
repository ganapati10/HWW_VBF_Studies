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
import warnings
#from inspect import signature
from funcsigs import signature
import bayesian_blocks


#                                                                                                                                                                                                          
# Sergio Blanco - 02/03/2020                                                                                                                                                                              
#                                                                                                                                                                                                           
# Adaptation of Bayesian Blocks from https://jakevdp.github.io/blog/2012/09/12/dynamic-programming-in-python/ for rebin a ROOT histogram.
#
# 
#


#                                                                                                                                                                                                          
# ------------------ Make new Bins for histograms ---------------------------------------                                                                                                                  
#                                                                                                                                                                                                         



class BayesianBlocks:

    _logger = logging.getLogger('BayesianBlocks')


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

        print_pop = True

    def makeBins(self, inputFile, outputDirPlots, variables, cuts, samples, plot, legend, groupplot):
            
        print "=================="
        print "==== makeBins ===="
        print "=================="

        
        
        ### All THE CONTENT OF THE CODE
        
        self._variables = variables
        self._samples   = samples
        self._cuts      = cuts
        
        ROOT.TH1.SetDefaultSumw2(True)
        
        ROOT.gROOT.cd()
        
        fileIn = ROOT.TFile(inputFile, "READ")
        
        for cutName in self._cuts:
            print "cut_", cutName
            
            if cutName!='VBF':
                continue

            for variableName, variable in self._variables.iteritems():
                
                if (variableName == 'mll_top'): continue

                print "========================="
                print "variable  ", variableName
                print "========================="
                print "Old range: ", variable['range']
                
                low = variable['range'][1]
                high = variable['range'][2]
                size = variable['range'][0]
                step = (high-low)/size
                
                content = np.zeros(size)
                place = np.linspace(low, high, size)
                total = 0
                
                for sampleName, plotdef in plot.iteritems():
                    
                    if 'samples' in variable and sampleName not in variable['samples']:
                        continue
                    if sampleName not in self._samples:
                        continue
                    if plotdef['isData']==1:
                        continue
                    
                    shapeName = cutName+"/"+variableName+"/histo_" + sampleName
                    #Check .root file                                                                                                                                                                      
                    if type(fileIn) is dict:
                        histo = fileIn[sampleName].Get(shapeName)     #Get the TH1 for each variable, cut and sample.                                                                                      
                    else:
                        histo = fileIn.Get(shapeName)
                    
                    histogram = histo.Clone("new_histo_" + sampleName + "_" + cutName + "_" + variableName)  #Open the .root file and create histogram
                    
                    for i in range(size):
                        content[i] = content[i] + histogram.GetBinContent(i)

                print "-----"
                print "The content of the Histogram is:  "
                print "-----"
                print content                    
            
                content = content + 1
                content = 2*np.round(content) #Factor two because of the low statistics //// Just for test: DON'T USE
                new_bins = bayesian_blocks.bayesian_blocks(place, content)
                
                print "========================="
                print "========================="
                print "-----Bayesian Blocks-----"
                print "New Binning for ", variableName
                print "========================="
                print "========================="
                
                print new_bins
                
                print ""
                print ""
                print "-------------------------------------------------------------------"
                


#                                                                                                                                                                                                           
#                                                                                                                                                                                                           
#   BIN MAKER                                                                                                                                                                                              
#                                                                                                                                                                                                           
#   INITIALIZATE THE BIN CLASS ABOVE                                                                                                                                                                       
#                                                                                                                                                                                                           
#                                                                                                                                                                                                           


if __name__=='__main__':
    sys.argv = argv

    print '''                                                                                                                                                                                               
    --------------------------------------------------------------------------------------------------                                                                                                      
        \  |         |                                                                                                                                                               
       |\/ |   _` |  |  /   _ \   __|                                                                                                                                                
       |   |  (   |    <    __/  |                                                                                                                                                   
      _|  _| \__,_| _|\_\ \___| _|                                                                                                                                                   
                                                                                                                                                                                                            
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

    if not opt.debug:
        pass
    elif opt.debug == 2:
        print 'Logging level set to DEBUG (%d)' % opt.debug
        logging.basicConfig( level=logging.DEBUG )
    elif opt.debug == 1:
        print 'Logging level set to INFO (%d)' % opt.debug
        logging.basicConfig( level=logging.INFO )


    plotter = BayesianBlocks()

    samples = OrderedDict()
    if opt.samplesFile == None:
        print " Please provide the samples structure (not strictly needed in mkPlot, since list of samples read from plot.py) "
    elif os.path.exists(opt.samplesFile):
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

    plotter.makeBins( opt.inputFile ,opt.outputDirPlots, variables, cuts, samples, plot, legend, groupPlot)

    print '.... Now clossing....'  
      
      
      
