from collections import OrderedDict
import common_configs

pt_edges = [200, 300, 400, 500, 650]
print(len(pt_edges))
configs = {
    "year": "UL17",
    "ModelName": "TTBarUL17",
    "gridHistFileName": "../Histograms/grid_oneScale.root",
    "histLocation": "../python/flat_templates/templates_UL17_1d.root",
    "binning": [50, 300, 5],
    "massScaleFactor": 2.0,
    "VaryOnlySignal": "True",
    "separateMassScales": "False",
    "pt_edges": pt_edges,
    "channels": OrderedDict(common_configs.top_channels(pt_edges)),
}
