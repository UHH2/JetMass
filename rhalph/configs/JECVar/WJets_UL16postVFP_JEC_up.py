from collections import OrderedDict
import common_configs

w_pt_edges = [500, 650, 800, 1200]
configs = {
    "year": "UL16postVFP",
    "ModelName": "WJetsUL16postVFPJECUP",
    "gridHistFileName": "../Histograms/grid_oneScale.root",
    "histLocation": "../python/flat_templates/templates_UL16postVFP_1d_jec_up.root",
    "binning": [50, 300, 5],
    "massScaleFactor": 2.0,
    "VaryOnlySignal": "True",
    "QCDSigmaScale": 10.0,
    "QCDFailConstant": "False",
    "JECVar":False,
    "separateMassScales": "False",
    "pt_edges": map(float, w_pt_edges),
    "channels": OrderedDict(common_configs.w_channels(w_pt_edges)),
}
