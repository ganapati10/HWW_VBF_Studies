combineCards.py ./datacards_2016_CP/hww2l2v_13TeV_SRVBF_HS/KD_H0M/datacard.txt > ./cards/H0M_HVV.txt
combineCards.py ./datacards_2016_CP/hww2l2v_13TeV_SRVBF_HS/KD_H0PH/datacard.txt > ./cards/H0PH_HVV.txt
combineCards.py ./datacards_2016_CP/hww2l2v_13TeV_SRVBF_HS/KD_H0L1/datacard.txt > ./cards/H0L1_HVV.txt


text2workspace.py cards/H0M_HVV.txt -o cards/H0M_HVV.root -P HiggsAnalysis.CombinedLimit.HWWCouplings:HWWCouplings --PO H0M > cards/scale.txt
text2workspace.py cards/H0PH_HVV.txt -o cards/H0PH_HVV.root -P HiggsAnalysis.CombinedLimit.HWWCouplings:HWWCouplings --PO H0PH > cards/scale.txt
text2workspace.py cards/H0L1_HVV.txt -o cards/H0L1_HVV.root -P HiggsAnalysis.CombinedLimit.HWWCouplings:HWWCouplings --PO H0L1 > cards/scale.txt


combine -d cards/H0M_HVV.root -n H0M_HVV -M MultiDimFit -t -1 -m 125 --algo grid --points 1000 --setParameters muV=1.0,Fai=0.0 --redefineSignalPOIs=Fai
combine -d cards/H0PH_HVV.root -n H0PH_HVV -M MultiDimFit -t -1 -m 125 --algo grid --points 1000 --setParameters muV=1.0,Fai=0.0 --redefineSignalPOIs=Fai
combine -d cards/H0L1_HVV.root -n H0L1_HVV -M MultiDimFit -t -1 -m 125 --algo grid --points 1000 --setParameters muV=1.0,Fai=0.0 --redefineSignalPOIs=Fai

rm combine_logger.out 

mv higgsCombineH0M_HVV.MultiDimFit.mH125*.root hists/higgsCombineH0M_HVV.MultiDimFit.mH125.root
mv higgsCombineH0PH_HVV.MultiDimFit.mH125*.root hists/higgsCombineH0PH_HVV.MultiDimFit.mH125.root
mv higgsCombineH0L1_HVV.MultiDimFit.mH125*.root hists/higgsCombineH0L1_HVV.MultiDimFit.mH125.root

python plotScan.py






