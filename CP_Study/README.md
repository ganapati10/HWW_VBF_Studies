# CP and HVV anomalous couplings for VBF Higgs production

This files are rewritten to run a likelihood scan over the fai parameters. To do that a new template must be done.

Reweighted the AC signals:

```
python ./Tools/TestSignalRW.py -b -l
```

Make plots from the new samples reweighted:

```
mkPlot.py --pycfg=configuration_CP.py --inputFile rootFileJJH/plots_JJH.root --showIntegralLegend 1
```


Make new templates for the scan:

```
python ./Tools/makeTemplates.py -b -l

mkDatacards.py --pycfg=configuration_hvv.py --inputFile rootFileJJH/plots_JJH_HVV.root
```

Then, move datacards to your combine folder and, for H0M, H0PH or H0l1, perform the scan:

```
combineCards.py hww2l2v_13TeV_SRVBF=./datacards_2016_CP/hww2l2v_13TeV_SRVBF/KD_H0M/datacard.txt > ./cards/H0M_HVV.txt

text2workspace.py cards/H0M_HVV.txt -o cards/H0M_HVV.root -P HiggsAnalysis.CombinedLimit.HWWCouplings:HWWCouplings --PO H0M > cards/scale.txt

combine -d cards/H0M_HVV.root -n H0M_HVV -M MultiDimFit -t -1 -m 125 --algo grid --points 1000 --setParameters muV=1.0,Fai=0.0 --redefineSignalPOIs=Fai

rm combine_logger.out

mv higgsCombineH0M_HVV.MultiDimFit.mH125*.root hists/higgsCombineH0M_HVV.MultiDimFit.mH125.root

python scripts/PlotScan.py
```





