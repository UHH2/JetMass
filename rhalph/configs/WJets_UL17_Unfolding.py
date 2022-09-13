from collections import OrderedDict
pt_edges = [500,575,650,725,800,1000,1200]
print(len(pt_edges))
configs  = {
  "ModelName":"VJetsUL17Unfolding",
  "gridHistFileName":"../Histograms/grid_oneScale.root",
  "histLocation":"../python/templates_UL17_1d_unfolding.root",
  "binning":[50,300,5],
  "massScaleFactor":2.0,
  "BernsteinOrders":[2,3],
  "InitialQCDFit":"False",
  "InitialQCDFitOrders":[1,3],
  "VaryOnlySignal":"False",
  "QCDSigmaScale":10.0,
  "QCDFailConstant":"False",
  "separateMassScales":"False",
  "pt_edges":[500.,575.,650.,725.,800.,1000.,1200.],
  "unfolding_bins":{
      "ptgen":[550, 650, 800, 1200,"Inf"],
      "msdgen":[0,55,65,72.5,80,87.5,95,110,125,150,180,220,"Inf"]
  },
  "channels":OrderedDict(
      [
          (
              "WJetsPt{}".format(pt_edges[ipt]),
              {
                  "QcdEstimation":"True",
                  "selection":"W",
                  "pt_bin":"{}to{}".format(*pt_edges[ipt:ipt+2]),
                  "samples":["QCD","WJetsMatched","WJetsUnmatched","ZJets","TTToHadronic","TTToSemiLeptonic"],
                  "NormUnc":{
                      # "WJets":1.2,
                      "ZJets":1.2
                  },
                  "signal":["WJetsMatched"],
                  "obs":"Data",
                  "regions":["pass","fail"]
              }
          )
          for ipt in range(len(pt_edges)-1)
  ])

}
