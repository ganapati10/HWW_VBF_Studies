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
        
        for variableName, variable in self._variables.iteritems():
            
            print "========================="
            print "variable  ", variableName
            print "========================="
            print "Old range: ", variable['range']
            
            low = variable['range'][1]
            high = variable['range'][2]
            size = variable['range'][0]
            step = (high-low)/size
            
            content = np.empty(size)
            place = np.linspace(low, high, size)
            
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
                    content[i] = histogram.GetBinContent(i)
            
            data = []
            j = 0
            for rep in content:
                for k in range(rep):
                    data.append(place[j])
                j = j + 1
            
            new_bins = bayesian_blocks(data)
            
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
      
      
      
      
      
      
      
    def bayesian_blocks(t):
        """Bayesian Blocks Implementation

        By Jake Vanderplas.  License: BSD
        Based on algorithm outlined in http://adsabs.harvard.edu/abs/2012arXiv1207.5578S

        Parameters
        ----------
        t : ndarray, length N
            data to be histogrammed

        Returns
        -------
        bins : ndarray
            array containing the (N+1) bin edges

        Notes
        -----
        This is an incomplete implementation: it may fail for some
        datasets.  Alternate fitness functions and prior forms can
        be found in the paper listed above.
        """
        # copy and sort the array
        t = np.sort(t)
        N = t.size

        # create length-(N + 1) array of cell edges
        edges = np.concatenate([t[:1],
                                0.5 * (t[1:] + t[:-1]),
                                t[-1:]])
        block_length = t[-1] - edges

        # arrays needed for the iteration
        nn_vec = np.ones(N)
        best = np.zeros(N, dtype=float)
        last = np.zeros(N, dtype=int)

        #-----------------------------------------------------------------
        # Start with first data cell; add one cell at each iteration
        #-----------------------------------------------------------------
        for K in range(N):
            # Compute the width and count of the final bin for all possible
            # locations of the K^th changepoint
            width = block_length[:K + 1] - block_length[K + 1]
            count_vec = np.cumsum(nn_vec[:K + 1][::-1])[::-1]

            # evaluate fitness function for these possibilities
            fit_vec = count_vec * (np.log(count_vec) - np.log(width))
            fit_vec -= 1.15  # 4 comes from the prior on the number of changepoints
            fit_vec[1:] += best[:K]

            # find the max of the fitness: this is the K^th changepoint
            i_max = np.argmax(fit_vec)
            last[K] = i_max
            best[K] = fit_vec[i_max]

        #-----------------------------------------------------------------
        # Recover changepoints by iteratively peeling off the last block
        #-----------------------------------------------------------------
        change_points =  np.zeros(N, dtype=int)
        i_cp = N
        ind = N
        while True:
            i_cp -= 1
            change_points[i_cp] = ind
            if ind == 0:
                break
            ind = last[ind - 1]
        change_points = change_points[i_cp:]

        return edges[change_points]            
     
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
      
      
      

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      




