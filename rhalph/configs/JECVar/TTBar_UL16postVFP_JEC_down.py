from collections import OrderedDict
import common_configs

top_pt_edges = [200, 300, 400, 500, 650]
configs = {
    "year": "UL16postVFP",
    "ModelName": "TTBarUL16postVFPJECDOWN",
    "gridHistFileName": "../Histograms/grid_oneScale.root",
    "histLocation": "../python/flat_templates/templates_UL16postVFP_1d_jec_down.root",
    "binning": [50, 300, 5],
    "massScaleFactor": 2.0,
    "BernsteinOrders": [2, 3],
    "InitialQCDFit": "False",
    "InitialQCDFitOrders": [2, 2],
    "VaryOnlySignal": "True",
    "QCDSigmaScale": 10.0,
    "QCDFailConstant": "False",
    "JECVar":False,
    "separateMassScales": "False",
    "pt_edges": map(float, top_pt_edges),
    "channels": OrderedDict(common_configs.top_channels(top_pt_edges)),
}
