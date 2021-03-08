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

class mkPlot2D:

    _logger = logging.getLogger('mkPlot2D')


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
        print "==== makePlot2D ===="
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
        
          #Loop under variables                                                                                                                                                                               
          for variableName, variable in self._variables.iteritems():

            if (len(variable['name']) == 1):
                continue
            
            tcanvas = ROOT.TCanvas("c_2D_" + variableName, "c_2D_" + variableName, 800, 600) #Iniciate the canvas, one for each variable.                                   
            
            thsData       = ROOT.THStack ("thsData_" + cutName + "_" + variableName,      "thsData_" + cutName + "_" + variableName)
            #print 'really before thstack ... one'
            thsSignal     = ROOT.THStack ("thsSignal_" + cutName + "_" + variableName,    "thsSignal_" + cutName + "_" + variableName)
            #print 'really before thstack ... two'
            thsBackground = ROOT.THStack ("thsBackground_" + cutName + "_" + variableName,"thsBackground_" + cutName + "_" + variableName)
            #print 'really before thstack ... three'

            
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
                    
                    histos[sampleName].SetMarkerColor(plotdef['color'])
                    histos[sampleName].SetMarkerSize(1)
                    histos[sampleName].SetMarkerStyle(20)
                    histos[sampleName].SetLineColor(self._getColor(plotdef['color']))
                    
                    # blind data
                    if 'isBlind' in plotdef.keys() and plotdef['isBlind'] == 1:
                      histos[sampleName].Reset()
                
                    # Per variable blinding
                    if 'blind' in variable:
                      if cutName in variable['blind']:
                        blind_range = variable['blind'][cutName]
                        if blind_range == "full":
                          for iBin in range(1, histos[sampleName].GetNbinsX()+1):
                            histos[sampleName].SetBinContent(iBin, 0)
                            histos[sampleName].SetBinError  (iBin, 0)
                          histos[sampleName].Reset()
                        elif type(blind_range) in [list,tuple] and len(blind_range)==2:
                          b0 = histos[sampleName].FindBin(blind_range[0])
                          b1 = histos[sampleName].FindBin(blind_range[1])
                          for iBin in range(1, histos[sampleName].GetNbinsX()+1):
                            if iBin >= b0 and iBin <= b1:
                              histos[sampleName].SetBinContent(iBin, 0)
                              histos[sampleName].SetBinError  (iBin, 0)
                    
                
                    thsData.Add(histos[sampleName])
                  
                else:
                  # MC style
                  # only background "filled" histogram
                  if plotdef['isSignal'] == 0:
                    #histos[sampleName].SetFillStyle(1001)
                    #histos[sampleName].SetFillColorAlpha(self._getColor(plotdef['color'],) 0.5)
                    #histos[sampleName].SetFillColor(self._getColor(plotdef['color']))
                    #histos[sampleName].SetLineColor(self._getColor(plotdef['color']+)1)
                    histos[sampleName].SetFillColor(self._getColor(plotdef['color']))
                    if 'fill' in plotdef:
                      histos[sampleName].SetFillStype(plotdef['fill'])
                    else:
                      histos[sampleName].SetFillStyle(3001)
                  else :
                    histos[sampleName].SetFillStyle(0)
                    histos[sampleName].SetLineWidth(2)
              
                  histos[sampleName].SetLineColor(self._getColor(plotdef['color']))
                  # scale to luminosity if MC
                  #histos[sampleName].Scale(self._lumi)  ---> NO! They are already scaled to luminosity in mkShape!
                  
                  if plotdef['isSignal'] == 1 :
                    if variable['divideByBinWidth'] == 1:
                      histos[sampleName].Scale(1,"width")
  
                    thsSignal.Add(histos[sampleName])
  
                  elif plotdef['isSignal'] == 2 or plotdef['isSignal'] == 3 :
                    #print "SigSup histo: ", histos[sampleName]
                    if  variable['divideByBinWidth'] == 1:
                      histos[sampleName].Scale(1,"width")
  
                    sigSupList.append(histos[sampleName])
  
                    if plotdef['isSignal'] == 3 :
                      #print "sigForAdditionalRatio histo: ", histos[sampleName]
                      sigForAdditionalRatioList[sampleName] = histos[sampleName]
                      sigForAdditionalDifferenceList[sampleName] = histos[sampleName]
                  else :
                    nexpected += histos[sampleName].Integral(1,histos[sampleName].GetNbinsX())   # it was (-1, -1) in the past, correct now
                    if variable['divideByBinWidth'] == 1:
                      histos[sampleName].Scale(1,"width")
  
                    thsBackground.Add(histos[sampleName])
                
                
          if len(groupPlot.keys()) == 0:          
            if thsBackground.GetNhists() != 0:
              thsBackground.Draw("hist same")
               
            if thsSignal.GetNhists() != 0:
              #for ihisto in range(thsSignal.GetNhists()) :
              #((thsSignal.GetHists().At(ihisto))).SetFillStyle(0)
              #((thsSignal.GetHists().At(ihisto))).Draw("hist same")
              thsSignal.Draw("hist same noclear")
          else :
            if thsBackground_grouped.GetNhists() != 0:
              thsBackground_grouped.Draw("hist same")
                 
            if thsSignal_grouped.GetNhists() != 0:
              thsSignal_grouped.Draw("hist same noclear")
              
            if len(sigSupList_grouped) != 0:
              for histo in sigSupList_grouped:
                histo.Draw("hist same")   
                
          
          
                    
                    
            
            
            
            
            
            
            
            
            
            
            
            
            
            

                #Just compute the variables equal to cuts: variable mjj with cut mjj_1                                                                                                                      
                if variableName != cutName.split('_')[0]:      #cutName must be variableName_i                                                                                                              
                    #print variableName + "not equal to " + cutName.split('_')[0]                                                                                                                           
                    continue

                print "variableName =", variableName

                # Define variables for store bin content of histograms                                                                                                                                      
                sig  = 0
                bkg  = 0
                data = 0

                print "Compute of Signal and Background processes"
                #print "Signal = 0"                                                                                                                                                                         
                #print "Background = 0"                                                                                                                                                                     

                ##-Loop above all samples in plot -> WW, Higgs, DY,..                                                                                                                                       
                for sampleName, plotdef in plot.iteritems():

                    if 'samples' in variable and sampleName not in variable['samples']:
                        continue

                    if sampleName not in self._samples:
                        continue


                    shapeName = cutName+"/"+variableName+"/histo_" + sampleName
                    print(shapeName)

                    #Check .root file                                                                                                                                                                       
                    if type(fileIn) is dict:
                        histo = fileIn[sampleName].Get(shapeName)     #Get the TH1 for each variable, cut and sample.                                                                                       
                    else:
                        histo = fileIn.Get(shapeName)

                    histogram = histo.Clone("new_histo_" + sampleName + "_" + cutName + "_" + variableName)  #Open the .root file and create histogram                                                      

                    if plotdef['isSignal'] == 1 :
                        sig = sig + histogram.GetEntries()

                    elif plotdef['isSignal'] == 0 :
                        bkg = bkg + histogram.GetEntries()

                    elif plotdef['isData'] == 1 :
                        data = data + histogram.GetEntries()

                    #print "------new sample " + sampleName + "------"                                                                                                                                      
                    #print "Signal = ", sig                                                                                                                                                                 
                    #print "Background = ", bkg                                                                                                                                                             

                ## End of samples loop.                                                                                                                                                                     

                if sig == 0 and bkg == 0:
                    print "Error: Non content in Signal/Background variables"
                    continue

                tHisto.SetBinContent(int(cutName.split("_")[1]), sig/math.sqrt(bkg+sig))   #Fill the histograms                                                                                             
                tpop.SetBinContent(int(cutName.split("_")[1]), sig)
                tpopB.SetBinContent(int(cutName.split("_")[1]), bkg)

                print "--------------------------------------------------------"
                print "New histogram computed, variable " + cutName
                print "Signal = ", sig
                print "Background = ", bkg
                print("Value of S/sqrt(S+B) : " + str(sig/math.sqrt(bkg+sig)))


            #End cuts loop                                                                                                                                                                                  

            #TH1F make up                                                                                                                                                                                   
            #tHisto.SetMinimum(0.0)                                                                                                                                                                         
            #tHisto.SetMaximum(1.0)                                                                                                                                                                         


            tHisto.SetMarkerStyle(ROOT.kFullCircle)
            tHisto.GetXaxis().SetTitle(variable['xaxis'])
            tHisto.GetYaxis().SetTitle('#S/\sqrt{B+S}')
            tHisto.SetTitle("")
            tHisto.Draw()

            legend = ROOT.TLegend(0.9, 0.87, 0.75, 0.82)

            legend.AddEntry(tHisto, variable['xaxis'])

            tHisto.SetStats(False)

            tHisto.Draw("p")

            #legend.Draw()                                                                                                                                                                                  

            tcanvas.Update()
            tcanvas.Draw()

            tcanvas2 = ROOT.TCanvas("Events " + variableName, "Events " + variableName, 800, 600)

            tpop.SetLineColor(ROOT.kRed)
            tpop.GetYaxis().SetTitle("Signal events per cut (Cumulative)")
            tpop.GetXaxis().SetTitle(variable['xaxis'])
            tpop.GetYaxis.SetLabelSize(0.02)
            tpop.SetStats(False)
            tpop.SetTitle("")
            tpop.Draw()
            tcanvas2.Update()

            rightmax = tpopB.GetMaximum()
            f1 = ROOT.TF1("f1","1", tpopB.GetXaxis().GetXmin(), tpopB.GetXaxis().GetXmax())
            if (rightmax!=0):
                scale = (tpop.GetMaximum() - tpop.GetMinimum())/rightmax
                tpopB.Scale(scale)
                tpopB.Add(f1, tpop.GetMinimum())

            tpopB.SetLineColor(ROOT.kBlue)
            tpopB.SetStats(False)
            tpopB.SetTitle("")
            tpopB.Draw("SAME")

            axis2 = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),0,rightmax,510,"+L")
            axis2.SetLineColor(ROOT.kBlue)
            axis2.SetTextColor(ROOT.kBlue)
            axis2.SetTitle("Background events per cut (Cumulative)")
            axis2.SetTitleOffset(1.)
            axis2.SetLabelSize(0.02)
            axis2.Draw()

            tcanvas2.Update()

            tcanvas2.Draw()

            tcanvas.SaveAs(self._outputDirPlots + "/" + "Sig_vs_Bkg_" + variableName + ".png")  #Save .png file                                                                                             
            tcanvas2.SaveAs(self._outputDirPlots + "/" + "Population_" + variableName + ".png")  #Save .png file                                                                                            
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


    plotter = mkPlot2D()

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
                
