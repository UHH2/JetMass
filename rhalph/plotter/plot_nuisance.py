import ROOT
from ROOT import gStyle,gROOT
gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
gStyle.SetOptTitle(0)

gStyle.SetTextFont(43)

gStyle.SetTitleOffset(0.86,"X")
gStyle.SetTitleOffset(1.8,"Y")
# gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadLeftMargin(0.19)
# gStyle.SetPadBottomMargin(0.15)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadTopMargin(0.08)
# gStyle.SetPadRightMargin(0.08)
gStyle.SetPadRightMargin(0.1)
gStyle.SetMarkerSize(0.5)
gStyle.SetHistLineWidth(2)
gStyle.SetTitleSize(0.02, "XYZ")
gStyle.SetLabelSize(0.04, "XYZ")
gStyle.SetNdivisions(506, "XYZ")
gStyle.SetLegendBorderSize(0)
gStyle.SetPadTickY(1)
gStyle.SetPadTickX(1)

def update_nuisance(name,fitargs,nuisances):
    nuisance = fitargs.find(name)
    effect = nuisances[name]['effect']
    nuisances[name]['central'] = 1 + nuisance.getVal() * effect
    nuisances[name]['down'] = 1 + (nuisance.getVal() - (nuisance.getAsymErrorLo() if nuisance.getAsymErrorLo() != 0 else nuisance.getErrorLo()) ) * effect
    nuisances[name]['up'] = 1 + (nuisance.getVal() + (nuisance.getAsymErrorHi() if nuisance.getAsymErrorHi() != 0 else nuisance.getErrorHi()) ) * effect
    
def plot_nuisances(config):

    f_r = ROOT.TFile(config['ModelName']+'/fitDiagnostics.root','READ')
    def dummy(effect=0.01):
        return {'effect':effect,'central':1,'up':1.1,'down':0.9}
        
    nuisances = {'wqq_norm':dummy(0.00089),'zqq_norm':dummy(0.00072),'CMS_lumi':dummy(0.027)}

    fitargs = f_r.Get('fit_s').floatParsFinal()

    for nu in nuisances.keys():
        update_nuisance(nu,fitargs,nuisances)
    print(nuisances)

def print_nuisance(config):
    f_r = ROOT.TFile(config['ModelName']+'/fitDiagnostics.root','READ')
    print(dir(f_r.Get('fit_s')))
    print(type(f_r.Get('fit_s').covQual()))
    c = ROOT.TCanvas("correlationHist","correlationHist",1000,1000)

    f_r.Get('fit_s').correlationHist().Draw("colz2")

    c.SaveAs('correlationHist.png')
    c.SaveAs('correlationHist.pdf')
    # fitargs = f_r.Get('fit_s').floatParsFinal()
    # for i in fitargs:
    #     print(i)
    # else:
    #     print('no more fitresults')

def plot_mass_scale_nuisances(config):
    f_grid = ROOT.TFile(config['gridHistFileName'])
    grid = f_grid.Get("grid")
    h_cat = f_grid.Get("categories");
    categories = [h_cat.GetXaxis().GetBinLabel(i) for i in range(1,h_cat.GetXaxis().GetNbins()+1)]
    f_r = ROOT.TFile(config['ModelName']+'/fitDiagnostics.root','READ')
    fitargs = f_r.Get('fit_s').floatParsFinal()


    Npt = grid.GetXaxis().GetNbins()
    Neta = grid.GetYaxis().GetNbins()
    Ncat = len(categories)
    variation = 0.1
    parname = []
    par = []
    errU = []
    errD = []
    for i in range(Npt):
        for j in range(Neta):
            for k in range(Ncat):
                parameter_name = "massScale_"
                parameter_name += "pt" + str(i)
                parameter_name += "_eta" + str(j)
                parameter_name += "_" + categories[k]
                params = fitargs.find(parameter_name)
                central = params.getValV()
                error_up = abs(params.getErrorHi())
                error_down = abs(params.getErrorLo())
                bincont = 1.0 + central * variation
                bincontU = 1.0 + (central+error_up) * variation
                bincontD = 1.0 + (central-error_up) * variation

                parname.append(parameter_name)
                par.append(central)
                errU.append(error_up)
                errD.append(error_down)

    Npoints = len(parname)
    g = ROOT.TGraphAsymmErrors(Npoints)
    for i in range(len(parname)):
        g.SetPoint(i, i+0.5, par[i])
        g.SetPointError(i, 0.0, 0.0, errD[i], errU[i])

    frame = ROOT.TH1F("frame", " ", len(parname), 0, len(parname))
    frame.GetYaxis().SetRangeUser(-2, 2)
    for i in range(len(parname)):
        frame.GetXaxis().SetBinLabel(i+1, parname[i])

    width = 200 + Npoints*75
    c = ROOT.TCanvas("c", "c", 2000, 600)
    ROOT.gPad.SetBottomMargin(.2)
    #gPad.SetRightMargin(.2)
    frame.GetYaxis().SetTitle("#sigma")
    frame.SetFillColor(ROOT.kWhite)
    frame.Draw()
    g.SetMarkerColor(ROOT.kBlack)
    g.SetMarkerStyle(8)
    g.SetMarkerSize(1)
    g.Draw("P SAME")
    c.SaveAs(config['ModelName']+'/plots/JetMassScale_Nuisance.pdf')



    
if __name__ == '__main__':
    import json
    config={'ModelName':'/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/rhalph/WMassModel'}
    #plot_nuisances(config)
    print_nuisance(config)
