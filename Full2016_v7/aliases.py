import os
import copy
import inspect

configurations = os.path.realpath(inspect.getfile(inspect.currentframe())) # this file
configurations = os.path.dirname(configurations) # ggH2016
configurations = os.path.dirname(configurations) # Differential
configurations = os.path.dirname(configurations) # Configurations
#configurations = '/afs/cern.ch/user/p/piedra/work/latinos/CMSSW_10_2_15_patch2/src/PlotsConfigurations/Configurations' # compute_SF.C

#aliases = {}

# imported from samples.py:
# samples, signals

mc = [skey for skey in samples if skey not in ('Fake', 'DATA', 'Dyemb')]
mc_emb = [skey for skey in samples if skey not in ('Fake', 'DATA')]

eleWP = 'mva_90p_Iso2016'
muWP = 'cut_Tight80x'

newEleWP = 'mva_90p_Iso2016' # same as nominal
#newEleWP = 'mva_90p_Iso2016_tthmva_70'
newMuWP = 'cut_Tight80x_tthmva_80'

aliases['LepWPCut'] = {
    ### eleWP + muWP combination
    'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP,
    ### eleWP + newMuWP combination
    #'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP+'*((abs(Lepton_pdgId[0])==11 || Muon_mvaTTH[Lepton_muonIdx[0]]>0.8) && (abs(Lepton_pdgId[1])==11 || Muon_mvaTTH[Lepton_muonIdx[1]]>0.8))',
    ### newEleWP + newMuWP combination
    #'expr': 'LepCut2l__ele_'+eleWP+'__mu_'+muWP+'*((abs(Lepton_pdgId[0])==11 || Muon_mvaTTH[Lepton_muonIdx[0]]>0.8) && (abs(Lepton_pdgId[1])==11 || Muon_mvaTTH[Lepton_muonIdx[1]]>0.8) && (abs(Lepton_pdgId[0])==13 || Electron_mvaTTH[Lepton_electronIdx[0]]>0.70) && (abs(Lepton_pdgId[1])==13 || Electron_mvaTTH[Lepton_electronIdx[1]]>0.70))',
    'samples': mc_emb + ['DATA']
}

aliases['gstarLow'] = {
    'expr': 'Gen_ZGstar_mass >0 && Gen_ZGstar_mass < 4',
    'samples': 'VgS'
}

aliases['gstarHigh'] = {
    'expr': 'Gen_ZGstar_mass <0 || Gen_ZGstar_mass > 4',
    'samples': 'VgS'
}

aliases['embedtotal'] = {
    #'expr': 'embed_total_WP90V1', For v7? Not in embed samples v6    wrt. eleWP
    'expr': 'embed_total_mva16',
    'samples': 'Dyemb'
}

# Fake leptons transfer factor
aliases['fakeW'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP,
    'samples': ['Fake']
}

# And variations - already divided by central values in formulas !
aliases['fakeWEleUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_EleUp',
    'samples': ['Fake']
}
aliases['fakeWEleDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_EleDown',
    'samples': ['Fake']
}
aliases['fakeWMuUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_MuUp',
    'samples': ['Fake']
}
aliases['fakeWMuDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_MuDown',
    'samples': ['Fake']
}
aliases['fakeWStatEleUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statEleUp',
    'samples': ['Fake']
}
aliases['fakeWStatEleDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statEleDown',
    'samples': ['Fake']
}
aliases['fakeWStatMuUp'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statMuUp',
    'samples': ['Fake']
}
aliases['fakeWStatMuDown'] = {
    'expr': 'fakeW2l_ele_'+eleWP+'_mu_'+muWP+'_statMuDown',
    'samples': ['Fake']
}

# gen-matching to prompt only (GenLepMatch2l matches to *any* gen lepton)
aliases['PromptGenLepMatch2l'] = {
    'expr': 'Alt$(Lepton_promptgenmatched[0]*Lepton_promptgenmatched[1], 0)',
    'samples': mc
}

aliases['Top_pTrw'] = {
    'expr': '(topGenPt * antitopGenPt > 0.) * (TMath::Sqrt(TMath::Exp(0.0615 - 0.0005 * topGenPt) * TMath::Exp(0.0615 - 0.0005 * antitopGenPt))) + (topGenPt * antitopGenPt <= 0.)',
    'samples': ['top']
}

# Jet bins
# using Alt$(CleanJet_pt[n], 0) instead of Sum$(CleanJet_pt >= 30) because jet pt ordering is not strictly followed in JES-varied samples

aliases['nCleanGenJet'] = {
    'linesToAdd': ['.L %s/Differential/ngenjet.cc+' % configurations],
    'class': 'CountGenJet',
    'samples': mc
}

# No jet with pt > 30 GeV
aliases['zeroJet'] = {
    'expr': 'Alt$(CleanJet_pt[0], 0) < 30.'
}

aliases['oneJet'] = {
    'expr': 'Alt$(CleanJet_pt[0], 0) > 30.'
}

aliases['multiJet'] = {
    'expr': 'Alt$(CleanJet_pt[1], 0) > 30.'
}

aliases['multiLepton'] = {
    'expr': 'Alt$(Lepton_pt[1], 0) > 30.'
}


aliases['mlljj'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/mlljj.cc+' ],
    'class' : 'mlljj'
}

aliases['mjjj'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/mjjj.cc+' ],
    'class' : 'mjjj'
}

aliases['jetdis'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/GetJetDis.cc+'],
    'class' : 'GetJetDis'
}

aliases['alpha1'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/alpha1.cc+'],
    'class' : 'alpha1'
}

aliases['alpha2'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/alpha2.cc+'],
    'class' : 'alpha2'
}

aliases['alpha3'] = {
    'linesToAdd' : ['.L /afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_4/src/PlotsConfigurations/Configurations/WW/Full2016_v6/extended/alpha3.cc+'],
    'class' : 'alpha3'
}


mes = ['RecoLevel_me_VBF_hsm', 'RecoLevel_me_VBF_hm', 'RecoLevel_me_VBF_hp', 'RecoLevel_me_VBF_hl', 'RecoLevel_me_VBF_mixhm', 'RecoLevel_me_VBF_mixhp', 'RecoLevel_me_VBF_mixhl', 'RecoLevel_me_VBF_S_hsm', 'RecoLevel_me_VBF_S_hm', 'RecoLevel_me_VBF_S_hp', 'RecoLevel_me_VBF_S_hl', 'RecoLevel_me_VBF_S_mixhm', 'RecoLevel_me_VBF_S_mixhp', 'RecoLevel_me_VBF_S_mixhl', 'RecoLevel_me_VBF_TU_hsm', 'RecoLevel_me_VBF_TU_hm', 'RecoLevel_me_VBF_TU_hp', 'RecoLevel_me_VBF_TU_hl', 'RecoLevel_me_VBF_TU_mixhm', 'RecoLevel_me_VBF_TU_mixhp', 'RecoLevel_me_VBF_TU_mixhl', 'RecoLevel_me_QCD_hsm', 'RecoLevel_me_QCD_hm', 'RecoLevel_me_QCD_hp', 'RecoLevel_me_QCD_hl', 'RecoLevel_me_QCD_mixhm', 'RecoLevel_me_QCD_mixhp', 'RecoLevel_me_QCD_mixhl']

mes2 = ['RecoLevel_me_QCD_hsm', 'RecoLevel_me_QCD_hm', 'RecoLevel_me_QCD_hp', 'RecoLevel_me_QCD_hl', 'RecoLevel_me_QCD_mixhm', 'RecoLevel_me_QCD_mixhp', 'RecoLevel_me_QCD_mixhl']

for me in mes:
    aliases[me]={
    'linesToAdd': [
    #'gSystem->Load("%s/src/ZZMatrixElement/MELA/data/%s/libmcfm_707.so","", kTRUE);'%(os.getenv('CMSSW_BASE'), os.getenv('SCRAM_ARCH')),
    #'gSystem->Load("libZZMatrixElementMELA.so","", kTRUE);',
    'gSystem->Load("%s/src/JHUGenMELA/MELA/data/%s/libmcfm_707.so","", kTRUE);'%(os.getenv('CMSSW_BASE'), os.getenv('SCRAM_ARCH')),
    'gSystem->Load("libJHUGenMELAMELA.so","", kTRUE);',
    '.L %s/patches/RecoLevelME_patch.cc+' % configurations],
    'class': 'RecoLevelME',
    'args': (me,)
    }


aliases['kd_smvbf']     = { 'expr': '1/(1+((RecoLevel_me_QCD_hsm)/RecoLevel_me_VBF_hsm))' }
aliases['kd_hmvbf']     = { 'expr': '1/(1+((RecoLevel_me_QCD_hsm)/(RecoLevel_me_VBF_hm)))' }
aliases['kd_hpvbf']     = { 'expr': '1/(1+((RecoLevel_me_QCD_hsm)/(RecoLevel_me_VBF_hp)))' }
aliases['kd_hlvbf']     = { 'expr': '1/(1+((RecoLevel_me_QCD_hsm)/(RecoLevel_me_VBF_hl)))' }
aliases['kd_vbf']       = { 'expr': 'max(max(kd_smvbf, kd_hmvbf), max(kd_hpvbf, kd_hlvbf))' }


##### DY Z pT reweighting
aliases['getGenZpt_OTF'] = {
    'linesToAdd':['.L %s/src/PlotsConfigurations/Configurations/patches/getGenZpt.cc+' % os.getenv('CMSSW_BASE')],
    'class': 'getGenZpt',
    'samples': ['DY']
}
handle = open('%s/src/PlotsConfigurations/Configurations/patches/DYrew30.py' % os.getenv('CMSSW_BASE'),'r')
exec(handle)
handle.close()
aliases['DY_NLO_pTllrw'] = {
    'expr': '('+DYrew['2016']['NLO'].replace('x', 'getGenZpt_OTF')+')*(nCleanGenJet == 0)+1.0*(nCleanGenJet > 0)',
    'samples': ['DY']
}
aliases['DY_LO_pTllrw'] = {
    'expr': '('+DYrew['2016']['LO'].replace('x', 'getGenZpt_OTF')+')*(nCleanGenJet == 0)+1.0*(nCleanGenJet > 0)',
    'samples': ['DY']
}

# B tagging

bAlgo = 'DeepB'
bWP = '0.2217'

aliases['bVeto'] = {
'expr': '(Sum$(CleanJet_pt > 20. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0) && mth > 60' }


aliases['bVetoDY'] = {                                                                                                              
'expr': '(Sum$(CleanJet_pt > 30. && abs(CleanJet_eta) < 2.5 && Jet_btagDeepB[CleanJet_jetIdx] > 0.2217) == 0) && mth < 60' }


aliases['bReq'] = {
'expr': '(   Alt$(CleanJet_pt[0],0) > 30. \
          && Alt$(CleanJet_pt[1],0) > 30. \
          && ( ( Alt$(abs(CleanJet_eta[0]),99)<2.5 && Alt$(Jet_btagDeepB[CleanJet_jetIdx[0]],0) > 0.2217 ) \
            || ( Alt$(abs(CleanJet_eta[1]),99)<2.5 && Alt$(Jet_btagDeepB[CleanJet_jetIdx[1]],0) > 0.2217 ) ) \
         )' 
}


### We need to define a CR for WW 

##aliases['wwcr'] = {
##    'expr': 'mth>60 && mtw2>30 && mll>100 && bVeto'
##}



# B tag scale factors

btagSFSource = '%s/src/PhysicsTools/NanoAODTools/data/btagSF/DeepCSV_2016LegacySF_V1.csv' % os.getenv('CMSSW_BASE')


aliases['Jet_btagSF_shapeFix'] = {
    'linesToAdd': [
        'gSystem->Load("libCondFormatsBTauObjects.so");',
        'gSystem->Load("libCondToolsBTau.so");',
        'gSystem->AddIncludePath("-I%s/src");' % os.getenv('CMSSW_RELEASE_BASE'),
        '.L %s/patches/btagsfpatch.cc+' % configurations
    ],
    'class': 'BtagSF',
    'args': (btagSFSource,),
    'samples': mc
}

aliases['bVetoSF'] = {
'expr': '( TMath::Exp(Sum$( TMath::Log( (CleanJet_pt>20 && abs(CleanJet_eta)<2.5)*Jet_btagSF_shapeFix[CleanJet_jetIdx]+1*(CleanJet_pt<20 || abs(CleanJet_eta)>2.5) ) ) ) )',
'samples': mc
}


aliases['bVetoDYSF'] = {                                                                                                            
'expr': '( TMath::Exp(Sum$( TMath::Log( (CleanJet_pt>30 && abs(CleanJet_eta)<2.5)*Jet_btagSF_shapeFix[CleanJet_jetIdx]+1*(CleanJet_pt<30 || abs(CleanJet_eta)>2.5) ) ) ) )',
'samples': mc
}                                                                                                                                    
aliases['bReqSF'] = {
'expr': '( ( ( Alt$(CleanJet_pt[0], 0)>30 && Alt$(abs(CleanJet_eta[0]),99)<2.5 )*( Alt$(Jet_btagSF_shapeFix[CleanJet_jetIdx[0]], 1) ) + ( Alt$(CleanJet_pt[0], 0)<30 || Alt$(abs(CleanJet_eta[0]),99)>2.5 ) )* \
           ( ( Alt$(CleanJet_pt[1], 0)>30 && Alt$(abs(CleanJet_eta[1]),99)<2.5 )*( Alt$(Jet_btagSF_shapeFix[CleanJet_jetIdx[1]], 1) ) + ( Alt$(CleanJet_pt[1], 0)<30 || Alt$(abs(CleanJet_eta[1]),99)>2.5 ) ) )\
        ',
'samples': mc
}


aliases['btagSF'] = {
    'expr': '( bVetoSF*bVeto +  bVetoDYSF*bVetoDY + bReqSF*bReq  + ( (!bVeto) && (!bVetoDY) &&  (!bReq) ) )',
'samples': mc
}



for shift in ['jes','lf','hf','lfstats1','lfstats2','hfstats1','hfstats2','cferr1','cferr2']:
    aliases['Jet_btagSF_shapeFix_up_%s' % shift] = {
        'class': 'BtagSF',
        'args': (btagSFSource, 'up_' + shift),
        'samples': mc
    }
    aliases['Jet_btagSF_shapeFix_down_%s' % shift] = {
        'class': 'BtagSF',
        'args': (btagSFSource, 'down_' + shift),
        'samples': mc
    }
    for targ in ['bVeto', 'bVetoDY', 'bReq']:
        alias = aliases['%sSF%sup' % (targ, shift)] = copy.deepcopy(aliases['%sSF' % targ])
        alias['expr'] = alias['expr'].replace('btagSF_shapeFix', 'btagSF_shapeFix_up_%s' % shift)

        alias = aliases['%sSF%sdown' % (targ, shift)] = copy.deepcopy(aliases['%sSF' % targ])
        alias['expr'] = alias['expr'].replace('btagSF_shapeFix', 'btagSF_shapeFix_down_%s' % shift)

    aliases['btagSF%sup' % shift] = {
        'expr': aliases['btagSF']['expr'].replace('SF', 'SF' + shift + 'up'),
        'samples': mc
    }

    aliases['btagSF%sdown' % shift] = {
        'expr': aliases['btagSF']['expr'].replace('SF', 'SF' + shift + 'down'),
        'samples': mc
    }

# data/MC scale factors
#aliases['new_SF'] = {   'linesToAdd': ['.L %s/patches/compute_SF.C+' % configurations],
#                        'class': 'compute_SF',
#                        'args' : ('2016', 2, 'total_SF'),
#                        'samples': mc
#}


#puidSFSource = '%s/src/LatinoAnalysis/NanoGardener/python/data/JetPUID_effcyandSF.root' % os.getenv('CMSSW_BASE')
puidSFSource = '{}/patches/PUID_80XTraining_EffSFandUncties.root'.format(configurations)

aliases['PUJetIdSF'] = {
    'linesToAdd': [
        'gSystem->AddIncludePath("-I%s/src");' % os.getenv('CMSSW_BASE'),
        '.L %s/patches/pujetidsf_event_new.cc+' % configurations
    ],
    'class': 'PUJetIdEventSF',
    'args': (puidSFSource, '2016', 'loose'),
    'samples': mc
}

aliases['ttHMVA_SF_2l'] = {'linesToAdd': ['.L %s/patches/compute_SF.C+' % configurations],
                           'class': 'compute_SF',
                           'args' : ('2016', 2, 'total_SF'),
                           'samples': mc_emb
}

# data/MC scale factors
aliases['SFweight'] = {
    'expr': ' * '.join(['SFweight2l', 'ttHMVA_SF_2l', 'LepWPCut', 'btagSF', 'PrefireWeight','PUJetIdSF']),
    #'expr': ' * '.join(['SFweight2l', 'LepSF2l__ele_' + eleWP + '__mu_' + muWP, 'LepWPCut', 'btagSF', 'PrefireWeight']),
    #'expr': ' * '.join(['SFweight2l', 'new_SF', 'LepWPCut', 'btagSF', 'PrefireWeight']),
    'samples': mc
}

# Muon ttHMVA SF needed for tau embedded samples
#aliases['Muon_ttHMVA_SF'] = {
#    'expr': '( (abs(Lepton_pdgId[0]) == 13)*(Lepton_tightMuon_cut_Tight80x_tthmva_80_IdIsoSF[0]/Lepton_tightMuon_cut_Tight80x_IdIsoSF[0])+(abs(Lepton_pdgId[0]) == 11) )*( (abs(Lepton_pdgId[1]) == 13)*(Lepton_tightMuon_cut_Tight80x_tthmva_80_IdIsoSF[1]/Lepton_tightMuon_cut_Tight80x_IdIsoSF[1])+ (abs(Lepton_pdgId[1]) == 11) )',
#    'samples' : ['Dyemb']
#}

# variations
aliases['SFweightEleUp'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__Up',
    'samples': mc_emb
}
aliases['SFweightEleDown'] = {
    'expr': 'LepSF2l__ele_'+eleWP+'__Do',
    'samples': mc_emb
}
aliases['SFweightMuUp'] = {
    'expr': 'LepSF2l__mu_'+muWP+'__Up',
    'samples': mc_emb
}
aliases['SFweightMuDown'] = {
    'expr': 'LepSF2l__mu_'+muWP+'__Do',
    'samples': mc_emb
}

aliases['nllWOTF'] = {
    'linesToAdd': ['.L %s/Differential/nllW.cc+' % configurations],
    'class': 'WWNLLW',
    'args': ('central',),
    'samples': ['WW']
}

# In WpWmJJ_EWK events, partons [0] and [1] are always the decay products of the first W
aliases['lhe_mW1'] = {
    'expr': 'TMath::Sqrt(2. * LHEPart_pt[0] * LHEPart_pt[1] * (TMath::CosH(LHEPart_eta[0] - LHEPart_eta[1]) - TMath::Cos(LHEPart_phi[0] - LHEPart_phi[1])))',
    'samples': ['WWewk']
}

# and [2] [3] are the second W
aliases['lhe_mW2'] = {
    'expr': 'TMath::Sqrt(2. * LHEPart_pt[2] * LHEPart_pt[3] * (TMath::CosH(LHEPart_eta[2] - LHEPart_eta[3]) - TMath::Cos(LHEPart_phi[2] - LHEPart_phi[3])))',
    'samples': ['WWewk']
}


# use HTXS_njets30 when moving to NanoAODv5 for all trees
#aliases['nCleanGenJet'] = {
#    'linesToAdd': ['.L %s/Differential/ngenjet.cc+' % configurations],
#    'class': 'CountGenJet',
#    'samples': signals
#}
#
## GGHUncertaintyProducer wasn't run for 2016 nAODv5 non-private
#thus = [
#    'ggH_mu',
#    'ggH_res',
#    'ggH_mig01',
#    'ggH_mig12',
#    'ggH_VBF2j',
#    'ggH_VBF3j',
#    'ggH_pT60',
#    'ggH_pT120',
#    'ggH_qmtop'
#]
#
#for thu in thus:
#    aliases[thu] = {
#        'linesToAdd': ['.L %s/Differential/gghuncertainty.cc+' % configurations],
#        'class': 'GGHUncertainty',
#        'args': (thu,),
#        'samples': ['ggH_hww'],
#        'nominalOnly': True
#    }


#thusQQ = [
#  "qqH_YIELD",
#  "qqH_PTH200",
#  "qqH_Mjj60",
#  "qqH_Mjj120",
#  "qqH_Mjj350",
#  "qqH_Mjj700",
#  "qqH_Mjj1000",
 # "qqH_Mjj1500",
 # "qqH_PTH25",
 # "qqH_JET01",
 # "qqH_EWK",
#]

#for thu in thusQQ:
#    aliases[thu] = {
#        'linesToAdd': ['.L %s/patches/qqhuncertainty.cc+' % configurations],
#        'class': 'QQHUncertainty',
#        'args': (thu,),
#        'samples': ['qqH_hww'],
#        'nominalOnly': True
#    }
    
