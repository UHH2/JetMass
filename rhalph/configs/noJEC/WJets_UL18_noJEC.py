from collections import OrderedDict
import common_configs

pt_edges = [500, 650, 800, 1200]
configs = {
    "year": "UL18",
    "ModelName": "WJetsUL18",
    "gridHistFileName": "../Histograms/grid_oneScale.root",
    "histLocation": "../python/flat_templates/templates_UL18_1d_jecpt.root",
    "binning": [50, 300, 5],
    "massScaleFactor": 2.0,
    "VaryOnlySignal": "True",
    "QCDSigmaScale": 10.0,
    "QCDFailConstant": "False",
    "separateMassScales": "False",
    "pt_edges": [500.0, 650.0, 800.0, 1200.0],
    "channels": OrderedDict(common_configs.w_channels(pt_edges)),
}
