from __future__ import print_function
import FitSubmitter
import json,argparse

if __name__ == "__main__":
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
