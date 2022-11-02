from collections import OrderedDict

pt_edges = [500, 575, 650, 725, 800, 1000, 1200]
print(len(pt_edges))
configs = {
    "ModelName": "VJetsUL17Unfolding",
    "gridHistFileName": "../Histograms/grid_oneScale.root",
    "histLocation": "../python/templates_UL17_1d_unfolding.root",
    "binning": [30, 300, 5],
    "massScaleFactor": 2.0,
    "BernsteinOrders": [2, 3],
    "InitialQCDFit": "False",
    "InitialQCDFitOrders": [1, 3],
    "VaryOnlySignal": "False",
    "QCDSigmaScale": 10.0,
    "QCDFailConstant": "False",
    "separateMassScales": "False",
    "pt_edges": [500.0, 575.0, 650.0, 725.0, 800.0, 1000.0, 1200.0],
    "unfolding_bins": {  # including overflow and underflow bin! those will be treated as background
        "ptgen": [0, 650, 800, 1200, "Inf"],
        "msdgen": [0.0, 44.5, 68.0, 80.5, 92.0, 132.5, "Inf"],
    },
    "channels": OrderedDict(
        [
            (
                "WJetsPt{}".format(pt_edges[ipt]),
                {
                    "QcdEstimation": "True",
                    "selection": "W",
                    "pt_bin": "{}to{}".format(*pt_edges[ipt: ipt + 2]),
                    "samples": ["QCD", "WJetsMatched", "WJetsUnmatched", "ZJets", "TTToHadronic", "TTToSemiLeptonic"],
                    "NormUnc": {
                        # "WJets":1.2,
                        "ZJets": 1.2
                    },
                    "signal": ["WJetsMatched"],
                    "obs": "Data",
                    "regions": ["pass", "fail"],
                },
            )
            for ipt in range(len(pt_edges) - 1)
        ]
    ),
}
