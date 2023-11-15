

def bernstein_orders(year, TF="Data", tagger="substructure", starting_at_30=False):
    if tagger == "":
        tagger = "substructure"
    # orders when fitting in [50,300]
    orders_50 = {
        "particlenet": {
            "Data": {
                # old orders before trigger fix
                # "UL16preVFP": [2, 5],
                # "UL16postVFP": [2, 5],
                # "UL17": [2, 6],
                # "UL18": [2, 6],
                "UL16preVFP": [2, 5],
                "UL16postVFP": [2, 5],
                "UL17": [2, 6],
                "UL18": [2, 6],
            },
            "QCD": {
                "UL16preVFP": [0, 2],
                "UL16postVFP": [1, 1],
                "UL17": [0, 1],
                "UL18": [0, 1],
            },
            "Data2TF": {
                "UL16preVFP": [2, 4],
                "UL16postVFP": [2, 3],
                "UL17": [2, 0],
                "UL18": [1, 0],
            },
        },
        "particlenetDDT": {
            "Data": {
                "UL16preVFP": [2, 3],
                "UL16postVFP": [2, 3],
                "UL17": [2, 3],
                "UL18": [2, 3],
            },
            "QCD": {},
            "Data2TF": {},
        },
        "substructure": {
            "Data": {
                # old orders before trigger fix
                # "UL16preVFP": [1, 4],
                # "UL16postVFP": [1, 2],
                # "UL17": [2, 3],
                # "UL18": [1, 4],
                "UL16preVFP": [0, 2],
                "UL16postVFP": [1, 2],
                "UL17": [2, 2],
                "UL18": [1, 4],
            },
            "QCD": {},
            "Data2TF": {},
        },
    }
    # orders when fitting in [30, 300]
    orders_30 = {
        "particlenet": {
            "Data": {
                "UL16preVFP": [2, 5],
                "UL16postVFP": [2, 5],
                "UL17": [2, 6],
                "UL18": [2, 6],
            },
        },
        "particlenetDDT": {
            "Data": {
                "UL16preVFP": [2, 3],
                "UL16postVFP": [2, 5],
                "UL17": [2, 1],
                "UL18": [2, 1],
            },
            "QCD": {},
            "Data2TF": {},
        },
        "substructure": {
            "Data": {
                "UL16preVFP": [1, 2],
                "UL16postVFP": [2, 3],
                "UL17": [2, 4],
                "UL18": [1, 5],
            },
            "QCD": {},
            "Data2TF": {},
        },
    }

    orders = orders_30 if starting_at_30 else orders_50
    return orders[tagger][TF][year]


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
