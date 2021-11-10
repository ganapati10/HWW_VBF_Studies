-- Load the library containing the matrix element
load_modules('/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_ME/build/libme_WW_ME.so')

wminus = declare_input("wminus")
wplus = declare_input("wplus")
jet1 = declare_input("jet1")
jet2 = declare_input("jet2")


parameters = {
    energy = 13000.,
}

inputs = {
    wminus.reco_p4,
    wplus.reco_p4,
    jet1.reco_p4,
    jet2.reco_p4,
    }

-- Build the partonic initial state
BuildInitialState.boost = {
    do_transverse_boost = true,
    particles = inputs,
}

-- Call matrix element on fully defined partonic event
MatrixElement.WW = {
    pdf = "CT10nlo",
    pdf_scale = 80.419002,
    matrix_element = "WW_ME_sm_P1_Sigma_sm_gg_wpwmuux",
    matrix_element_parameters = {
        card = "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/WW_ME/Cards/param_card.dat"
    },
    initialState = "boost::partons",
    particles = {
        inputs = inputs,
        ids = {
            {
                pdg_id = 24,
                me_index = 1,
            },
            {
                pdg_id = -24,
                me_index = 2,
            },
            {
                pdg_id = 1,
                me_index = 3,
            },
	    {
                pdg_id = -1,
                me_index = 4,
            },
        }
    },
}

-- Define quantity to be returned to MoMEMta
integrand("WW::output")











