#!/usr/bin/env python
import ROOT
import os
from CombineUtils import CombineWorkspace as cw


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-w", "--workspace", type=str, help="path to workspace ROOT-file. Should contain workspace names 'w'."
    )
    parser.add_argument("-t", "--toysFile", type=str, help="path to ROOT-file containing toys.")
    parser.add_argument("-o", "--outdir", type=str, default="toyPlots/")
    args = parser.parse_args()

    ws = cw(args.workspace)

    ws.load_toys(args.toysFile)

    ROOT.gROOT.SetBatch(True)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    for cat in ws.cats:
        c = ROOT.TCanvas(cat)
        frame = ws.obs.frame()
        ws.data_sets[cat].plotOn(frame)
        for iToy, toy in enumerate(ws.toys):

            toy_color = 46 + iToy
            toy[cat].plotOn(
                frame,
                ROOT.RooFit.LineColor(toy_color),
                ROOT.RooFit.MarkerColor(toy_color),
                ROOT.RooFit.LineStyle(2),
                ROOT.RooFit.MarkerStyle(26),
            )
        frame.Draw()

        outdir = args.outdir + "/"#toy%i/"%iToy
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        c.SaveAs(outdir + cat + ".pdf")
        del c, frame
