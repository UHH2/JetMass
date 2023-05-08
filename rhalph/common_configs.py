

def bernstein_orders(year, particlenet=False, TF="Data"):
    orders_particlenet = {
        "Data": {
            "UL16preVFP": [2, 5],
            "UL16postVFP": [2, 4],
            "UL17": [2, 4],
            "UL18": [2, 4],
        },
        "QCD": {},
        "Data2TF": {},
    }
    orders_substructure = {
        "Data": {
            "UL16preVFP": [2, 2],
            "UL16postVFP": [1, 2],
            "UL17": [2, 3],
            "UL18": [2, 3],
        },
        "QCD": {},
        "Data2TF": {},
    }
    if particlenet:
        return orders_particlenet[TF][year]
    else:
        return orders_substructure[TF][year]


def w_channels(w_pt_edges):
    return [
        (
            "WJetsPt{}".format(w_pt_edges[ipt]),
            {
                "QcdEstimation": "True",
                "selection": "W",
                "pt_bin": "{}to{}".format(*w_pt_edges[ipt:ipt + 2]),
                "samples": ["QCD", "WJetsMatched", "WJetsUnmatched", "ZJets", "TTToHadronic", "TTToSemiLeptonic"],
                "NormUnc": {
                    # "WJets": {"value": 1.2, "decorrelateRegions": True},
                    "WJets": 1.2,
                    "ZJets": 1.2
                },
                "signal": ["WJetsMatched"],
                "obs": "Data",
                "regions": ["pass", "fail"],
            },
        )
        for ipt in range(len(w_pt_edges) - 1)
    ]


def top_channels(top_pt_edges):
    return [
        (
            "TopPt{}".format(top_pt_edges[ipt]),
            {
                "selection": "top",
                "pt_bin": "{}to{}".format(*top_pt_edges[ipt: ipt + 2]),
                "samples": [
                    "QCD",
                    "ST_s",
                    "ST_t",
                    "ST_tW",
                    "DYJets",
                    "WJets",
                    "TTToHadronic",
                    "TTTo2L2Nu",
                    "TTToSemiLeptonic_mergedTop",
                    "TTToSemiLeptonic_mergedW",
                    "TTToSemiLeptonic_mergedQB",
                    "TTToSemiLeptonic_notMerged",
                ],
                "NormUnc": {
                    "QCD": 2.0,
                    "DYJets": 2.0,
                    "ST": 1.23,
                    "WJets": 1.19,
                    "TTTo": {"value": 1.20, "decorrelateRegions": True},
                },
                "signal": ["TTToSemiLeptonic_mergedTop", "TTToSemiLeptonic_mergedW"],
                "obs": "Data",
                "regions": ["pass", "passW", "fail"],
            },
        )
        for ipt in range(len(top_pt_edges) - 1)
    ]
