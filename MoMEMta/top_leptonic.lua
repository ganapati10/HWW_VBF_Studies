-- Load the library containing the matrix element
load_modules('/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/top_leptonic/build/libme_top_leptonic.so')

wminus = declare_input("wminus")
wplus = declare_input("wplus")
bjet1 = declare_input("bjet1")
bjet2 = declare_input("bjet2")


parameters = {
    energy = 13000.,
}

inputs = {
    wminus.reco_p4,
    wplus.reco_p4,
    bjet1.reco_p4,
    bjet2.reco_p4,
    }

-- Build the partonic initial state
BuildInitialState.boost = {
    do_transverse_boost = true,
    particles = inputs,
}

-- Call matrix element on fully defined partonic event
MatrixElement.ttbar = {
    pdf = "CT10nlo",
    pdf_scale = 173.,
    matrix_element = "top_leptonic_sm_P1_Sigma_sm_gg_wpbwmbx",
    matrix_element_parameters = {
        card = "/afs/cern.ch/work/s/sblancof/public/CMSSW_10_6_10/top_leptonic/Cards/param_card.dat"
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
                pdg_id = 5,
                me_index = 2,
            },
            {
                pdg_id = -24,
                me_index = 3,
            },
	    {
                pdg_id = -5,
                me_index = 4,
            },
        }
    },
}

-- Define quantity to be returned to MoMEMta
integrand("ttbar::output")











