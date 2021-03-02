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
# Sergio Blanco - 02/03/2020                                                                                                                                                                                
#                                                                                                                                                                                                           


#                                                                                                                                                                                                           
# ------------------ Make new Bins for histograms ---------------------------------------                                                                                                                   
#                                                                                                                                                                                                           

class reBin:

    _logger = logging.getLogger('reBin')


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

    def makePlot(self, inputFile, outputDirPlots, variables, cuts, samples, plot, legend, groupplot):
            
        print "=================="
        print "==== makeBins ===="
        print "=================="


      
      ### All THE CONTENT OF THE CODE
      
      self._variables = variables
      self._samples   = samples
      self._cuts      = cuts
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
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


    plotter = Plot_Sig_Bkg()

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

    plotter.makePlot( opt.inputFile ,opt.outputDirPlots, variables, cuts, samples, plot, legend, groupPlot)

    print '.... Now clossing....'  
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      




