from __future__ import print_function
import FitSubmitter
import json,argparse

if __name__ == "__main__":
    import fitplotter
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=str, help="path to json with config")
    parser.add_argument('--justplots',action='store_true', help='just make plots.')
    parser.add_argument('--build', action='store_true', help='just build workspace for combine')
    parser.add_argument('--job_index', default="", type=str)
    parser.add_argument('--minimalModel',action='store_true')
    parser.add_argument('--noMassScales',action='store_true')
    parser.add_argument('--defaultPOI',action='store_true')
    parser.add_argument('--skipTemplatePlots',action='store_true')
    parser.add_argument('--customCombineWrapper',action='store_true')
    parser.add_argument('--combineOptions',type=str,help="string with additional cli-options passed to combine",default="")
    # parser.add_argument('--QCDOnly',action="store_true",help="perform full TF fit to QCD MC - for BernsteinOptimization only")
    args = parser.parse_args()


    submitter = FitSubmitter.FitSubmitter("WJets.json","DUST/test_dir",dry_run=False)
    submitter.display()
    print('1TF:')
    submitter.scan_TF_orders( (range(1,3),range(1,3)) , combine_method = "GOF")
    LumiScan = {'Pseudo':[
        ("L4fb",['lumiScale:0.1']),
        ("L41fb",['lumiScale:1.0']),
        ("L100fb",['lumiScale:2.3889154324']),
        ("L200fb",['lumiScale:4.7778308648'])]}

    submitter.scan(LumiScan,"LumiScan",combine_method = "GOF")
