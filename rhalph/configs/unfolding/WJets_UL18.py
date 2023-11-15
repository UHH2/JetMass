from collections import OrderedDict
import common_configs

# pt_edges = [500, 575, 650, 725, 800, 1000, 1200, "Inf"]
pt_edges = [575, 650, 725, 800, 1000, 1200, "Inf"]
configs = {
    "year": "UL18",
    "ModelName": "WJetsUL18Unfolding",
    "gridHistFileName": "../Histograms/grid_oneScale.root",
    "histLocation": "../python/flat_templates/templates_UL18_1d_unfolding.root",
    # "Pseudo": ["fromMC"],
    "binning": [50, 300, 5],
    "massScaleFactor": 2.0,
    "VaryOnlySignal": "False",
    "QCDSigmaScale": 10.0,
    "QCDFailConstant": "False",
    "separateMassScales": "False",
    "regularization": ["msd", "pt"],
    "regularizationStrength": 0.825,
    "uniformGenbins": "True",
    "scaleGenBinWidth": "False",
    # "pt_edges": [500.0, 575.0, 650.0, 725.0, 800.0, 1000.0, 1200.0, "Inf"],
    "pt_edges": [575.0, 650.0, 725.0, 800.0, 1000.0, 1200.0, "Inf"],
    "pt_cutoff": 1400,
    "msd_cutoff": 100,
    # "xsec_priors": {"WJetsMatched": 0.734689},
    # "xsec_priors": {"WJetsMatched": 0.702038},
    "unfolding_bins": {  # including overflow and underflow bin! those will be treated as background
        "ptgen": [500, 650, 800, 1200, "Inf"],
        "msdgen": [0., 70., 80., 90., "Inf"],
    },
    "channels": OrderedDict(common_configs.w_channels(pt_edges)),
}
