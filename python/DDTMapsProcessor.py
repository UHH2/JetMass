#!/usr/bin/env coffea_python.sh
import awkward as ak
from coffea.nanoevents import BaseSchema
from coffea import processor, hist
from coffea.analysis_tools import PackedSelection
from coffea.util import save
import os
import glob
from coffea_util import CoffeaWorkflow

import numpy as np


class DDTMapPrep(processor.ProcessorABC):
    def __init__(self):
        hists = {}

        year_axis = hist.Cat("year", "Year")

        pT_ax = hist.Bin("pt", "$p_{T}$ [GeV]", 500, 0.0, 5000.0)
        rho_ax = hist.Bin("rho", "$\\rho$ [GeV]", 300, -10.0, 0.0)
        N2_ax = hist.Bin("n2", "N2", 100, 0, 1.5)

        hists.update(
            {"leadingJet": hist.Hist("Events", year_axis, pT_ax, rho_ax, N2_ax)}
        )
        hists.update(
            {
                "leadingJet_corrected_pt_mass": hist.Hist(
                    "Events", year_axis, pT_ax, rho_ax, N2_ax
                )
            }
        )
        hists.update(
            {
                "leadingJet_corrected_pt": hist.Hist(
                    "Events", year_axis, pT_ax, rho_ax, N2_ax
                )
            }
        )

        self._hists = processor.dict_accumulator(
            {
                "nevents": processor.defaultdict_accumulator(float),
                "sumw": processor.defaultdict_accumulator(float),
                "sumw2": processor.defaultdict_accumulator(float),
                **hists,
            }
        )

    @property
    def accumulator(self):
        return self._hists

    def process(self, events):
        out = self.accumulator.identity()

        year = events.metadata["dataset"]

        selection = PackedSelection()
        selection.add("cleaner", (events.pt > 170) & (abs(events.eta) < 2.4))
        selection.add("jetpfid", events.jetpfid == 1)
        selection.add("trigger", events['trigger_bits'][:, 7] == 1)

        events = events[selection.require(cleaner=True, jetpfid=True, trigger=True)]

        jecfactor = events.jecfactor

        # pt = events.pt
        # pt_raw = pt/jecfactor
        pt_raw = events.pt
        pt = pt_raw * jecfactor

        mjet_raw = events.mjet
        mjet = mjet_raw * jecfactor

        rho = 2 * np.log(mjet_raw / pt_raw)
        out["leadingJet"].fill(
            year=year, pt=pt_raw, rho=rho, n2=events.N2, weight=events.weight
        )

        rho_corrected_pt = 2 * np.log(mjet_raw / pt)
        out["leadingJet_corrected_pt"].fill(
            year=year, pt=pt, rho=rho_corrected_pt, n2=events.N2, weight=events.weight
        )

        rho_corrected_pt_mass = 2 * np.log(mjet / pt)
        out["leadingJet_corrected_pt_mass"].fill(
            year=year,
            pt=pt,
            rho=rho_corrected_pt_mass,
            n2=events.N2,
            weight=events.weight,
        )

        out["nevents"][year] += len(events)
        out["sumw"][year] += ak.sum(events.weight)
        out["sumw2"][year] += ak.sum(events.weight * events.weight)

        return out

    def postprocess(self, accumulator):
        return accumulator


if __name__ == "__main__":
    workflow = CoffeaWorkflow("DDTMaps")

    workflow.parser.add_argument(
        "--output", "-o", type=str, default="qcd_pt_v_rho_v_n2.coffea"
    )
    workflow.parser.add_argument("--years", nargs="+", default=["UL17"])

    args = workflow.parse_args()

    workflow.processor_instance = DDTMapPrep()
    workflow.processor_schema = BaseSchema

    # path = '/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/vjetsTrees/workdir_{SELECTION}_{YEAR}/'
    path = "/nfs/dust/cms/user/albrechs/UHH2/JetMassOutput/vjetsTrees/ForDDTMaps/workdir_{SELECTION}_ddt_{YEAR}/"

    sample_pattern = os.path.join(path, "*QCD*.root")

    samples = {
        y: {"files": glob.glob(sample_pattern.format(SELECTION=args.selection, YEAR=y))}
        for y in args.years
    }

    for k, v in samples.items():
        print(k, len(v["files"]))

    output_file_path = os.path.join(os.getcwd(), args.output)

    os.chdir(os.environ["TMPDIR"])
    if args.scaleout > 0:
        workflow.init_dask_htcondor_client(1, 10, 5)

    output = workflow.run(samples)

    if not output_file_path.endswith(".coffea"):
        output_file_path += ".coffea"

    save(output, output_file_path)
