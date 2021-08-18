import sys
import ROOT 
import numpy as np
import shutil
import math
from os import path

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

#src = path.realpath("rootFileJJH/plots_JJH.root")
src = path.realpath("rootFile/plots_WW_2016_CP.root")
#########################################################


def CompareHWWRW(Cat, Var, Prod, Hyp):

 print " " 
 print Cat, Var, Prod, Hyp
 print " " 

 f = ROOT.TFile.Open(''+src+'', 'read')

 if Hyp == "H0PM" : H1 = f.Get('hww2l2v_13TeV_'+Cat+'/'+Var+'/histo_'+Prod+'H0PM') 
 else             : H1 = f.Get('hww2l2v_13TeV_'+Cat+'/'+Var+'/histo_'+Prod+'H0PM_'+Hyp+'') 
 if Hyp == "H0M" :  H2 = f.Get('hww2l2v_13TeV_'+Cat+'/'+Var+'/histo_'+Prod+'H0M') 
 else            :  H2 = f.Get('hww2l2v_13TeV_'+Cat+'/'+Var+'/histo_'+Prod+'H0M_'+Hyp+'') 
 if Hyp == "H0Mf05" :  H3 = f.Get('hww2l2v_13TeV_'+Cat+'/'+Var+'/histo_'+Prod+'H0Mf05') 
 else               :  H3 = f.Get('hww2l2v_13TeV_'+Cat+'/'+Var+'/histo_'+Prod+'H0Mf05_'+Hyp+'') 

 H1.SetDirectory(0)
 H2.SetDirectory(0)
 H3.SetDirectory(0)

 f.Close()

 ########################## Naive sum ###############################
 
 Sum = H1.Clone()
 Sum.SetDirectory(0)
 Sum.Add(H2,1)
 Sum.Add(H3,1)
 Sum.Scale(0.3333)

 ################# weighted sum using Average bit #####################

 SumA = H1.Clone()
 SumA.SetDirectory(0)

 H1A = H1.Clone()
 H2A = H2.Clone()
 H3A = H3.Clone()

 H1A.SetDirectory(0)
 H2A.SetDirectory(0)
 H3A.SetDirectory(0)

 H1A.SetBit(ROOT.TH1.kIsAverage)
 H2A.SetBit(ROOT.TH1.kIsAverage)
 H3A.SetBit(ROOT.TH1.kIsAverage)

 H1A.Add(H2A,1)
 H1A.Add(H3A,1)

 for i in range(1, H1A.GetXaxis().GetNbins()+1):
   n = H1A.GetBinContent(i)
   e = H1A.GetBinError(i)
   SumA.SetBinContent(i, n)
   SumA.SetBinError(i, e)

 #############################

 H1.SetLineColor(ROOT.kRed)
 H2.SetLineColor(ROOT.kGreen)
 H3.SetLineColor(ROOT.kOrange)

 Sum.SetLineColor(ROOT.kGreen)
 SumA.SetLineColor(ROOT.kRed)

 H1.SetLineWidth(1)
 H2.SetLineWidth(1)
 H3.SetLineWidth(1)

 Sum.SetLineWidth(3)
 SumA.SetLineWidth(3)

 canvasM = ROOT.TCanvas('canvasM', '', 500, 500)

 H1.SetMinimum(-2*Sum.GetMaximum())
 H1.SetMaximum(2*Sum.GetMaximum())
 H1.GetXaxis().SetTitle(""+Var+"")

 H1.Draw("e")
 H2.Draw("same e")
 H3.Draw("same e")
 Sum.Draw("same hist")
 SumA.Draw("same hist")

 legend = ROOT.TLegend(0.75,0.6,1.0,1.0)
 legend.AddEntry(H1,"H0PM","l")
 legend.AddEntry(H2,"H0M","l")
 legend.AddEntry(H3,"H0Mf05","l")
 legend.AddEntry(Sum,"Sum","l")
 legend.AddEntry(SumA,"Avg Sum","l")
 legend.Draw()

 # canvasM.SetLogy()
 canvasM.SaveAs("plotJJH/HGGRW_"+Cat+"_"+Var+"_"+Prod+Hyp+".pdf")
 canvasM.SaveAs("plotJJH/HGGRW_"+Cat+"_"+Var+"_"+Prod+Hyp+".png")

##########################################################


VBFJJConfig = [ ("SRVBF", "D_int_hm_2", "VBF_", "H0M"),
                ("SRVBF", "D_int_hp_2", "VBF_", "H0PH"),
                ("SRVBF", "D_int_hm_2", "VBF_", "H0L1"),          
]

SigConfig = VBFJJConfig
for cat, var, prod, hyp in SigConfig :
 CompareHWWRW(cat, var, prod, hyp)
