#!/nfs/dust/cms/user/albrechs/miniconda3/envs/coffea/bin/python

from coffea.util import load,save
def flatten_templates(hists,hist_name,samples,jec_applied_on='pt&mJ'):
    hist_axes = hists[hist_name].axes
    templates = {}
        
    for pt_bin in range(len(hist_axes['pt'].edges)):
        pt_inclusive = pt_bin==len(hist_axes['pt'].edges)-1
        if(pt_inclusive):
            pt_bin_name = 'inclusive'
        else:
            pt_bin_name = f"{hist_axes['pt'].edges[pt_bin]}to{hist_axes['pt'].edges[pt_bin+1]}"
            pt_bin_name = pt_bin_name.replace('.0','')
            pt_bin_name = pt_bin_name.replace('inf','Inf')
        for sample in samples:
            if(sample not in hist_axes['dataset']):
                continue
            if(sample == 'Data' and 'variation' in hist_name):
                continue
            legacy_hist_name = hist_name.replace('vjets_',f'W_{sample}__')
            legacy_hist_name = legacy_hist_name.replace('ttbar_',f'top_{sample}__')
            legacy_hist_name = legacy_hist_name.replace('variation','')
            legacy_hist_name = legacy_hist_name.replace('mjet_',f'mjet_{pt_bin_name}_')
            if(pt_inclusive):
                h_ = hists[hist_name][:,::sum,sample,jec_applied_on]
            else:
                h_ = hists[hist_name][:,pt_bin,sample,jec_applied_on]
            
            templates[legacy_hist_name] = h_
    return templates


if(__name__ == "__main__"):
    import argparse
    from utils import hist_to_th1
    import uproot

    parser = argparse.ArgumentParser()
    parser.add_argument("--input","-i",default="templates_2017.coffea")
    parser.add_argument("--output","-o",default="templates_1d_2017")
    parser.add_argument("--JEC",default = "pt&mJ", choices=['none','pt','pt&mJ'])
    
    args = parser.parse_args()
    
    hists = load(args.input)

    selections = {
        'vjets':{
            'regions':['inclusive','pass','fail'],
            'samples':['Data','WJetsMatched','WJetsUnmatched','ZJetsMatched','ZJetsUnmatched','QCD'],
        },
        'ttbar':{
            'regions':['inclusive','pass','passW','fail'],
            'samples':[],
        }
    }
    print(f'flattening templates and saving to ROOT from {args.input}')
    templates = {}
    for selection in selections.keys():
        samples = selections[selection]['samples']
        for region in selections[selection]['regions']:
            templates.update(flatten_templates(hists,f'{selection}_mjet_{region}',samples=samples, jec_applied_on=args.JEC))
            for variation in ['all','chargedH','neutralH','gamma','other']:
                templates.update(flatten_templates(hists,f'{selection}_mjet_variation_0_0_{variation}_{region}__up',samples=samples, jec_applied_on=args.JEC))
                templates.update(flatten_templates(hists,f'{selection}_mjet_variation_0_0_{variation}_{region}__down',samples=samples, jec_applied_on=args.JEC))
    save(templates,f"{args.output}.coffea")
    fout = uproot.recreate(f"{args.output}.root")
    for hist_name,H in templates.items():
        fout[hist_name] = hist_to_th1(H,hist_name)
    fout.close()
        
