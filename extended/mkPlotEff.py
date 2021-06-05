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
# ------------------ make 2D plots ---------------------------------------                                                                                                                   
#                                                                                                                                                                                                          

class mkPlotEff:

    _logger = logging.getLogger('mkPlotEff')


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
      
        print "===================="
        print "==== makePlotEff ===="
        print "===================="

        self._variables = variables
        self._samples   = samples
        self._cuts      = cuts

        self._outputDirPlots = outputDirPlots
        os.system ("mkdir " + outputDirPlots + "/")

        ROOT.TH2.SetDefaultSumw2(True)

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
        
        for cutName in self._cuts:
            print "cut_", cutName
            
            if (cutName != "DNN"): continue
            
            #Loop under variables                                                                                                                                                                          
            for variableName, variable in self._variables.iteritems():
                
                if (len(variable['name']) == 2):
                    continue
                 
                if (variable['name'] != 'DNN_vbf'): continue
                    
                tcanvas = ROOT.TCanvas("c_Eff_" + variableName, "c_Eff_" + variableName, 200, 200, 800, 600) #Iniciate the canvas, one for each variable.                                   
                
                rango = variable['range']

                thsData       = ROOT.THStack ("thsData_" + cutName + "_" + variableName,      "thsData_" + cutName + "_" + variableName)
                #print 'really before thstack ... one'
                thsSignal     = ROOT.THStack ("thsSignal_" + cutName + "_" + variableName,    "thsSignal_" + cutName + "_" + variableName)
                #print 'really before thstack ... two'
                thsBackground = ROOT.THStack ("thsBackground_" + cutName + "_" + variableName,"thsBackground_" + cutName + "_" + variableName)
                #print 'really before thstack ... three'

                histos = {}

                for sampleName, plotdef in plot.iteritems():
                    if 'samples' in variable and sampleName not in variable['samples']:
                        continue
                    shapeName = cutName+"/"+variableName+'/histo_' + sampleName
                    print '     -> shapeName = ', shapeName
                    if type(fileIn) is dict:
                        histo = fileIn[sampleName].Get(shapeName)
                    else:
                        histo = fileIn.Get(shapeName)
                    print ' --> ', histo
                    print 'new_histo_' + sampleName + '_' + cutName + '_' + variableName
                    histos[sampleName] = histo.Clone('new_histo_' + sampleName + '_' + cutName + '_' + variableName)

                    if plotdef['isData']==1:
                        continue
                        
                    histogram2D = ROOT.TH2F(sampleName+"_EffHist_"+variableName, sampleName+"_EffHist_"+variableName, len(np.asarray(rango[0]))-1, np.asarray(rango[0]), len(np.asarray(rango[1]))-1, np.asarray(rango[1]))
                    print "------------"
                    print "Bin test:"
                    for i in range(1, len(np.asarray(rango[0]))):
                        
                        for j in range(1, len(np.asarray(rango[1]))):
                            index = (len(np.asarray(rango[0]))-1)*(i-1)+j
                            histogram2D.SetBinContent(i, j, histos[sampleName].GetBinContent((len(np.asarray(rango[0]))-1)*(i-1)+j))
                            print "", index
                    print "-------------"
                            
                            
                    #histogram2D.SetFillColor(self._getColor(plotdef['color']))
                    #histogram2D.SetFillStyle(3001)
                    
                    if sampleName == "qqH_hww":
                        histogram2D.SetMinimum(0)
                        histogram2D.GetXaxis().SetTitle("pt_{ll}")
                        histogram2D.GetYaxis().SetTitle("m_{jj}")
                        histogram2D.SetTitle("qqH 2D Distribution")
                        histogram2D.SetStats(False)
                        histogram2D.Draw("COLZ")
                        tcanvas.SetLogz(False)
                        tcanvas.Draw()
                        tcanvas.SaveAs(self._outputDirPlots + "/" + cutName + "_c_2D_" + variableName + "_" + sampleName + ".png")  #Save .png file
                    
                    elif sampleName == "ggH_hww":
                        
                        histogram2D.SetMinimum(0)
                        histogram2D.GetXaxis().SetTitle("pt_{ll}")
                        histogram2D.GetYaxis().SetTitle("m_{jj}")
                        histogram2D.SetTitle("ggH 2D Distribution")
                        histogram2D.SetStats(False)
                        histogram2D.Draw("COLZ")
                        tcanvas.SetLogz(False)
                        tcanvas.Draw()
                        tcanvas.SaveAs(self._outputDirPlots + "/" + cutName + "_c_2D_" + variableName + "_" + sampleName + ".png")  #Save .png file
                    else:
                        continue
                    


    def _getColor(self, color):
      if type(color) == int:
        return color
      elif type(color) == tuple:
        # RGB
        return ROOT.TColor.GetColor(*color)
      elif type(color) == str:
        # hex string
        return ROOT.TColor.GetColor(color)



#                                                                                                                                                                                                          
#                                                                                                                                                                                                          
#   PLOT MAKER                                                                                                                                                                                             
#                                                                                                                                                                                                          
#   INITIALIZATE THE PLOT CLASS ABOVE                                                                                                                                                                      
#                                                                                                                                                                                                          
#                                                                                                                                                                                                          


if __name__=='__main__':
    sys.argv = argv

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

    if not opt.debug:
        pass
    elif opt.debug == 2:
        print 'Logging level set to DEBUG (%d)' % opt.debug
        logging.basicConfig( level=logging.DEBUG )
    elif opt.debug == 1:
        print 'Logging level set to INFO (%d)' % opt.debug
        logging.basicConfig( level=logging.INFO )


    plotter = mkPlotEff()

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
