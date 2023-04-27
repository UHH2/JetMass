import FitSubmitter
from CombineWorkflows import CombineWorkflows
import json

# dirty hack to keep this primarily python3 compatible, but still be able to run in py2 wihtout changes
try:
    range = xrange
except BaseException as e:
    print(e)
    pass


def submit_ftest_prep(
    config,
    ntoys,
    njobs,
    QCDOrders=None,
    DataOrders=None,
    algo="saturated",
    dry_run=False,
    PseudoData=False,
    unfolding=False,
):
    for i in range(njobs):
        do_qcd_test = False
        do_2TF = QCDOrders is not None

        if (
            QCDOrders is not None and (isinstance(QCDOrders[0], range) or isinstance(QCDOrders[0], list))
        ) or DataOrders is None:
            do_qcd_test = True
        do_data_test = False
        if (
            not do_qcd_test
            and DataOrders is not None
            and (isinstance(DataOrders[0], range) or isinstance(DataOrders[0], list))
        ):
            do_data_test = True
        if do_qcd_test and do_data_test:
            raise ValueError(
                "You provided fixed orders for both QCD and Data TF. You have to provide a lists of orders that should"
                " be scanned for either TF!"
            )
        if not do_qcd_test and not do_data_test:
            raise ValueError(
                "You provided running orders for both QCD and Data TF. You have to provide lists of orders that should"
                " be scanned for just one TF!"
            )

        if PseudoData or do_qcd_test:
            config["Pseudo"] = [""]
        else:
            config["Pseudo"] = []

        # config["QCDSigmaScale"] = 0.01

        cw = CombineWorkflows()
        cw.combineCMSSW = (
            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_6_28/CMSSW_10_6_28/src/UHH2/JetMass/rhalph/CMSSW_10_2_13/"
        )
        cw.method = "FTestBatch"
        cw.seed = 42 + i
        cw.POI = "r"
        cw.freezeParameters = "r"
        cw.extraOptions = (
            ("--toysFrequentist " if do_data_test else "")
            + ("--toysNoSystematics " if do_qcd_test else "")
            # ("--toysFrequentist --setParameters r=1 " if do_data_test else "")
            # + ("--toysNoSystematics --setParameters r=0 " if do_qcd_test else "")
            + " --cminDefaultMinimizerStrategy 2 --cminDefaultMinimizerTolerance 0.01"
        )
        # +(
        #     ""
        #     if do_qcd_test
        #     else " --redefineSignalPOIs massScale_pt0_eta0_all_PT500,massScale_pt0_eta0_all_PT650,"
        # "massScale_pt0_eta0_all_PT800"
        # )

        cw.algo = algo
        # cw.toysOptions = ("--toysFrequentist --expectSignal 1" if do_data_test else "") + (
        #     "--toysNoSystematics --expectSignal 0" if do_qcd_test else ""
        # )
        cw.toysOptions = ""
        cw.toys = ntoys

        workdir = (
            "FTest_"
            + config["ModelName"]
            + "_"
            + ("QCDTFScan" if do_qcd_test else "")
            + ("DataTFScan" + ("_QCDOrder%ix%i" % QCDOrders if do_2TF else "") if do_data_test else "")
            + "_%s_Seed%i" % (cw.algo, cw.seed)
        )

        fs = FitSubmitter.FitSubmitter(config, "DUST/" + workdir, dry_run=dry_run)
        fs.batchname = workdir.replace("_", "")
        fs.CombinePath = cw.combineCMSSW
        # fs.extra_jetmass_options = "--build --noMassScales"
        fs.extra_jetmass_options = "--build" + (" --unfolding --one-bin" if unfolding else "")
        fs.fit_qcd_model = do_qcd_test
        fs.scan_TF_orders(DataOrders, QCDOrders, combine_method=cw)


def ftest_QCDMC_TF_fitDiagnostics_closure(config, QCDOrders=(range(0, 7), range(0, 7)), dry_run=False):
    cw = CombineWorkflows()
    cw.set_method("diagnostics")
    cw.POI = "r"
    cw.freezeParameters = ""
    workdir = "FitDiagClosure_1TFScan"

    fs = FitSubmitter.FitSubmitter(config, "DUST/" + workdir, dry_run=dry_run)
    fs.batchname = workdir.replace("_", "")
    fs.extra_jetmass_options = "--build"  # just build workspace and not perform fit with main workspace
    fs.fit_qcd_model = True

    fs.do_fit_diagnostics_plots()

    fs.scan_TF_orders(([0], [0]), QCDOrders, combine_method=cw)


def ftest_Data_TF_fitDiagnostics_closure(
        config, Orders=(range(0, 7), range(0, 7)), QCDOrders=([3], [1]),
        PseudoData=False, dry_run=False
):
    if PseudoData:
        config["Pseudo"] = [""]
    else:
        config["Pseudo"] = []

    cw = CombineWorkflows()
    cw.set_method("diagnostics")
    cw.POI = "r"
    cw.freezeParameters = ""
    workdir = "FitDiagClosure_2TFScan_%sData_QCDOrders_%ix%i" % (
        "Pseudo" if PseudoData else "",
        QCDOrders[0][0],
        QCDOrders[1][0],
    )

    fs = FitSubmitter.FitSubmitter(config, workdir, dry_run=dry_run)
    fs.batchname = workdir.replace("_", "")
    fs.extra_jetmass_options = "--build --noMassScales"  # just build workspace and not perform fit with main workspace

    fs.do_fit_diagnostics_plots()

    fs.scan_TF_orders(Orders, QCDOrders, combine_method=cw)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="path to json file holding the fit-config")
    parser.add_argument(
        "--ntoys", type=int, default=50, help="specify how many toy events should be generated for each job."
    )
    parser.add_argument(
        "--njobs", type=int, default=10, help="Specify how many jobs should be submitted for each FTest config."
    )
    parser.add_argument("--dryrun", action="store_true", help="perfrom a dry-run")
    parser.add_argument("--unfolding", action="store_true", help="Ftests for unfolding setup. (onegebin)")
    args = parser.parse_args()

    if args.config.endswith(".json"):
        configs = json.load(open(args.config))
    else:
        execfile(args.config)

    # Orders = (range(0, 5)), (range(0, 5))
    Orders = (range(0, 7)), (range(0, 7))

    for algo in ["saturated"]:
        # for algo in ["saturated"]:
        # if algo == "saturated":
            # # FTests with QCD-TF orders 0x4 and scanning Data-TF orders
            # submit_ftest_prep(
            #     config, args.ntoys, args.njobs, QCDOrders=(0, 4), DataOrders=Orders, algo=algo, dry_run=args.dryrun
            # )
            # # FTests without QCD-TF and scanning Data-TF orders
        submit_ftest_prep(
            configs,
            args.ntoys,
            args.njobs,
            QCDOrders=None,
            DataOrders=Orders,
            algo=algo,
            dry_run=args.dryrun,
            unfolding=args.unfolding,
        )
        # FTests scanning QCD-TF orders
        # submit_ftest_prep(config, args.ntoys, args.njobs, QCDOrders=Orders,  algo=algo, dry_run=args.dryrun)
