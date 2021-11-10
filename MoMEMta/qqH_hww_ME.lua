-- Load the library containing the matrix element
load_modules('/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/qqH_hww_ME/build/libme_qqH_hww_ME.so')

higgs = declare_input("higgs")
jet1 = declare_input("jet1")
jet2 = declare_input("jet2")


parameters = {
    energy = 13000.,
}

inputs = {
    higgs.reco_p4,
    jet1.reco_p4,
    jet2.reco_p4,
    }

-- Build the partonic initial state
BuildInitialState.boost = {
    do_transverse_boost = true,
    particles = inputs,
}

-- Call matrix element on fully defined partonic event
MatrixElement.qqh = {
    pdf = "CT10nlo",
    pdf_scale = 125.,
    matrix_element = "qqH_hww_ME_HPO_ewk_prod_NLO_P1_Sigma_HPO_ewk_prod_NLO_uu_huu",
    matrix_element_parameters = {
        card = "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/qqH_hww_ME/Cards/param_card.dat"
    },
    initialState = "boost::partons",
    particles = {
        inputs = inputs,
        ids = {
            {
                pdg_id = 25,
                me_index = 1,
            },
            {
                pdg_id = 1,
                me_index = 2,
            },
            {
                pdg_id = -1,
                me_index = 3,
            },
        }
    },
}

-- Define quantity to be returned to MoMEMta
integrand("qqh::output")











