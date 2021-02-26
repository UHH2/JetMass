from __future__ import print_function
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
import glob,os,warnings
    
def nlls(filename):
    try:
        f = ROOT.TFile(filename)
        limits = f.Get('limit')
        nll=[l.limit for l in limits]
        status=[l.limitErr for l in limits]
        return nll,status
    except:
        # raise BaseException("One of the Fits seems to have failed!")
        #warnings.simplefilter("once", UserWarning)
        warnings.warn("Not able to retrieve limits from "+filename+" !")
        return [None],[None]

class FTest():
    def __init__(self):
        self._nbins=15*5
        self._maxNToys = -1
        self._orders_1 = (0,0)
        self._orders_2 = (0,0)
        
        self._nll_base1 = None
        self._nll_base2 = None

        self._nll_toys1 = None
        self._nll_toys2 = None

        self._F_base = None
        self._F_toys = None

        self._nToys = -1

        self._algo = "saturated"
        self._lumi = 41.8
        self._isData = False
        self._outDir = "."
        self.results = {}

        self._npall = 1
        
    @property
    def orders_1(self):
        return self._orders_1

    @orders_1.setter
    def orders_1(self,orders):
        self._orders_1 = orders
        self._np1 = (orders[0]+1)*(orders[1]+1) + self._npall

    @property
    def orders_2(self):
        return self._orders_2

    @orders_2.setter
    def orders_2(self,orders):
        self._orders_2 = orders
        self._np2 = (orders[0]+1)*(orders[1]+1) + self._npall

    def prob(self, nll_base, nll_toys):
        return float(sum(nll_base > nll_toys))/float(len(nll_toys))

    def fstat(self,nll1,nll2):
        result = ((nll1-nll2)/(self._np2-self._np1)) / (nll2/(self._nbins-self._np2))
        if('numpy.ndarray' in str(type(nll1))):
            result=result[np.where(nll1>nll2)]
        # else:
        #     if(nll1>nll2):
        #         print('problematic?\t',"nll1=%f,nll2=%f,np1=%f,np2=%f,nbins=%f"%(nll1,nll2,self._np1,self._np2,self._nbins))
        #         return None
        # commented this out since it would result into crashs as i did not implement another check further into ftest plotting/prob calc.
        # Will leave this as it is, but keep in mind how negative observed F-vals are possible (i.e. nll1<nll2)
        return result
        
    def process_files(self,workdir,dir_name_template,skip_failed_fits=False):
        """
        skip_failed_fits: if True fits where limitErr != 0 are skipped. 
        (per default limitErr is 0 when using GoodnessOfFit methods in combine. So to use this, one needs to tweak the used combine installation)
        """
        itoys_max = 1e9  # *100
        toy_workdirs = glob.glob(workdir)
                
        file_dir_1 = (dir_name_template%(self.orders_1[0],self.orders_1[1])).replace("/qcdmodel","")+'/'+dir_name_template%(self.orders_1[0],self.orders_1[1])
        file_dir_2 = (dir_name_template%(self.orders_1[0],self.orders_1[1])).replace("/qcdmodel","")+'/'+dir_name_template%(self.orders_2[0],self.orders_2[1])

        file_path_1 = ""
        file_path_2 = ""
        for toy_workdir in toy_workdirs:
            seed = toy_workdir.split("Seed")[-1].split('/')[0]
            file_path_1 = toy_workdir+'/'+file_dir_1+'/higgsCombineBaseline.GoodnessOfFit.mH0.%s.root'%seed
            file_path_2 = toy_workdir+'/'+file_dir_2+'/higgsCombineBaseline.GoodnessOfFit.mH0.%s.root'%seed
            if(os.path.isfile(file_path_1) and os.path.isfile(file_path_2)):
                self._nll_base1 = nlls(file_path_1)[0][0]
                self._nll_base2 = nlls(file_path_2)[0][0]
                if(self._nll_base1 is not  None and self._nll_base2 is not None):
                    break
        if(self._nll_base1 is  None or self._nll_base2 is None):
            warnings.warn("Could not find files for both fits. skipping (%i,%i).."%self.orders_1)
            print("Tried following paths last: \n %s \n %s"%(file_path_1,file_path_2))
            return 1

                

        if(self._nll_base1 is None or self._nll_base2 is None):
            warnings.warn("One of the Baseline fits failed. skipping")
            return 1
        
        nll_toys1 = []
        nll_toys2 = []
        status_toys1 = []
        status_toys2 = []
        for i,w in enumerate(toy_workdirs):
            if(i>itoys_max):
                break
            seed = w.split("Seed")[-1].split('/')[0]
            toy_file_1 = w+'/'+file_dir_1+'/higgsCombine.GoodnessOfFit.mH0.%s.root'%seed
            toy_file_2 = w+'/'+file_dir_2+'/higgsCombine.GoodnessOfFit.mH0.%s.root'%seed
            
            this_nll_toys1,this_status_toys1 = nlls(toy_file_1)
            this_nll_toys2,this_status_toys2 = nlls(toy_file_2)

            if((this_nll_toys1[0] is None) or (this_nll_toys2[0] is None) or (len(this_nll_toys1) != len(this_nll_toys2)) ):
                print("Toys from seed dir %s are empty or not compatible."%seed)
                continue
            
            nll_toys1 += this_nll_toys1
            nll_toys2 += this_nll_toys2
            status_toys1 += this_status_toys1
            status_toys2 += this_status_toys2
            
        if(skip_failed_fits):
            status_1 = np.array(status_toys1)
            status_2 = np.array(status_toys2)
            converged_fits = np.where((status_1>=0) & (status_2>=0))
            if(len(converged_fits[0])==0):
                print("All Fits faild..skipping")
                return -1
            self._nll_toys1 = np.array(nll_toys1)[converged_fits]
            self._nll_toys2 = np.array(nll_toys2)[converged_fits]
            print(len(nll_toys1),len(self._nll_toys1))
        else:
            self._nll_toys1 = np.array(nll_toys1)
            self._nll_toys2 = np.array(nll_toys2)
        if(len(nll_toys1) == len(nll_toys2)):
            self._nToys = len(nll_toys1)
        else:
            warnings.warn("Number of Toys does not match between compared TF-Order configurations. skipping..")
            return -1
            
        
        self._F_base = self.fstat(self._nll_base1, self._nll_base2)
        self._F_toys = self.fstat(self._nll_toys1, self._nll_toys2)

        if(len(self._F_toys) == 0):
            warnings.warn("F_vals for toys is emtpy (Probably all negative). skipping..")
            return -1
        return 0

        
    # def plot(self, nll_base, nll_toys, label=""):
    def plot_ftest(self):
        ftest_label = "FTest_TF%ix%i_v_TF%ix%i"%(*self.orders_1,*self.orders_2)
        ftest_prob = self.prob(self._F_base, self._F_toys)
        self.results.update({ftest_label:ftest_prob})
        print(ftest_label,"prob:",ftest_prob)        
        self.plot(self._F_toys, self._F_base, ftest_prob, ftest_label, "FTest")

    def plot_gofs(self):
        if(self._nll_toys1 is None or len(self._nll_toys1) == 0):
            print("nll1 empty. No GOF plot will be created!")
            return
        GOF1_label = "GOF_TF%ix%i"%self.orders_1
        GOF1_prob = self.prob(self._nll_base1, self._nll_toys1)
        print(GOF1_label,"prob:",GOF1_prob)
        self.plot(self._nll_toys1, self._nll_base1, GOF1_prob, GOF1_label, "GoodnessOfFit")

        # GOF2_label = "GOF_TF%ix%i"%self.orders_2
        # GOF2_prob = self.prob(self._nll_base2, self._nll_toys2)
        # print(GOF2_label,"prob:",GOF2_prob)
        # self.plot(self._nll_toys2, self._nll_base2, GOF2_prob, GOF2_label, "GoodnessOfFit")

    
    def plot(self,iToys,iCentral,prob,iLabel,method):
        """
        This plotting function is strongly (almost entirely) based on https://github.com/kakwok/ZPrimePlusJet/blob/50bb522fd9cac6f195f38bcb1c98c9f510168e4d/fitting/PbbJet/limit.py#L108
        """
        lCan   = ROOT.TCanvas(str(iLabel),str(iLabel),800,600)    
        lCan.SetLeftMargin(0.12) 
        lCan.SetBottomMargin(0.12)
        lCan.SetRightMargin(0.1)
        lCan.SetTopMargin(0.1)

        max_x = max(max(iToys),iCentral)
        if method=='FTest':
            lH = ROOT.TH1F(iLabel+"hist",iLabel,70,0,max_x+1)
            lH_cut = ROOT.TH1F(iLabel+"hist",iLabel.replace("_"," "),70,0,max_x+1)
        elif method=='GoodnessOfFit' and self._algo=='saturated':
            lH = ROOT.TH1F(iLabel+"hist",iLabel,70,0,max_x+100)
            lH_cut = ROOT.TH1F(iLabel+"hist",iLabel.replace("_"," "),70,0,max_x+100)
        elif method=='GoodnessOfFit' and self._algo=='KS':
            lH = ROOT.TH1F(iLabel+"hist",iLabel.replace("_"," "),70,0,max_x+0.05)
            lH_cut = ROOT.TH1F(iLabel+"hist",iLabel+"hist",70,0,max_x+0.05)
        elif method=='GoodnessOfFit' and self._algo=='AD':
            lH = ROOT.TH1F(iLabel+"hist",iLabel.replace("_"," "),70,0,max_x+10)
            lH_cut = ROOT.TH1F(iLabel+"hist",iLabel+"hist",70,0,max_x+10)
    
        if method=='FTest':
            lH.GetXaxis().SetTitle("F = #frac{-2log(#lambda_{1}/#lambda_{2})/(p_{2}-p_{1})}{-2log#lambda_{2}/(n-p_{2})}")
            lH.GetXaxis().SetTitleSize(0.025)
            lH.GetXaxis().SetTitleOffset(2)
            lH.GetYaxis().SetTitle("Pseudodatasets")
            lH.GetYaxis().SetTitleOffset(0.85)
        elif method=='GoodnessOfFit' and self._algo=='saturated':
            lH.GetXaxis().SetTitle("-2log#lambda")  
            lH.GetYaxis().SetTitle("Pseudodatasets")
            lH.GetYaxis().SetTitleOffset(0.85)
        elif method=='GoodnessOfFit' and self._algo=='KS':
            lH.GetXaxis().SetTitle("KS")  
            lH.GetYaxis().SetTitle("Pseudodatasets")
            lH.GetYaxis().SetTitleOffset(0.85)
        elif method=='GoodnessOfFit' and self._algo=='AD':
            lH.GetXaxis().SetTitle("AD")  
            lH.GetYaxis().SetTitle("Pseudodatasets")
            lH.GetYaxis().SetTitleOffset(0.85)
        for val in iToys:
            lH.Fill(val)
            if val > iCentral:
                lH_cut.Fill(val)
        lH.SetMarkerStyle(20)
        lH.Draw("pez")
        lLine  = ROOT.TArrow(iCentral,0.25*lH.GetMaximum(),iCentral,0)
        lLine.SetLineColor(ROOT.kBlue+1)
        lLine.SetLineWidth(2)
        
        lH_cut.SetLineColor(ROOT.kViolet-10)
        lH_cut.SetFillColor(ROOT.kViolet-10)
        lH_cut.Draw("histsame")
    
        if method=='FTest':
            fdist = ROOT.TF1("fDist", "[0]*TMath::FDist(x, [1], [2])", 0,max_x+1)
            fdist.SetParameter(0,lH.Integral()*((max_x+1)/70.))
            fdist.SetParameter(1,self._np2-self._np1)
            fdist.SetParameter(2,self._nbins-self._np2)
            fdist.Draw('same')
        elif method=='GoodnessOfFit' and self._algo=='saturated':
            chi2_func = ROOT.TF1('chisqpdf','[0]*ROOT::Math::chisquared_pdf(x,[1])',0,max_x+100)
            chi2_func.SetParameter(0,lH.Integral())
            chi2_func.SetParameter(1,50)
            chi2_func.Draw('same')
            lH.Fit(chi2_func,"mle")        
        lH.Draw("pezsame")
        lLine.Draw()
        
        tLeg = ROOT.TLegend(0.6,0.6,0.89,0.89)
        tLeg.SetLineColor(ROOT.kWhite)
        tLeg.SetLineWidth(0)
        tLeg.SetFillStyle(0)
        tLeg.SetTextFont(42)
        tLeg.AddEntry(lH,"toy data (Ntoys=%i)"%len(iToys),"lep")
        tLeg.AddEntry(lLine,"observed = %.3f"%iCentral,"l")
        tLeg.AddEntry(lH_cut,"p-value = %.3f"%(1-prob),"f")
        if method=='FTest':
            tLeg.AddEntry(fdist,"F-dist, ndf = (%.0f, %.0f) "%(fdist.GetParameter(1),fdist.GetParameter(2)),"l")        
        elif method=='GoodnessOfFit' and self._algo=='saturated':
            tLeg.AddEntry(chi2_func,"#chi^{2} fit, ndf = %.1f #pm %.1f"%(chi2_func.GetParameter(1),chi2_func.GetParError(1)),"l")
        
        tLeg.Draw("same")

        l = ROOT.TLatex()
        l.SetTextAlign(11)
        l.SetTextSize(0.06)
        l.SetTextFont(62)
        l.SetNDC()
        l.DrawLatex(0.12,0.91,"CMS")
        l.SetTextSize(0.05)
        l.SetTextFont(52)
        if self._isData:
            l.DrawLatex(0.23,0.91,"Preliminary")
        else:
            l.DrawLatex(0.23,0.91,"Simulation")
        l.SetTextFont(42)
        l.DrawLatex(0.70,0.91,"%.1f fb^{-1} (13 TeV)"%self._lumi)
        l.SetTextFont(52)
        l.SetTextSize(0.045)
    
        lCan.SaveAs(self._outDir+'/'+iLabel+".pdf")


    def print_results(self):
        for label,prob in self.results.items():
            print(label+' p-value:\t%.2f'%(1-prob))



def process_TF_prep(dir_pattern,model_dir_pattern,output_dir,algo='KS',n_p_other=1):
    ft = FTest()
    ft._algo = algo
    ft._outDir = output_dir
    ft._npall = n_p_other
    if(not os.path.isdir(ft._outDir)):
        os.makedirs(ft._outDir)
    for i_pt in range(0,6):
        for i_rho in range(0,6):
            for (j_pt,j_rho) in [(i_pt+1,i_rho),(i_pt,i_rho+1)]:
                print("#"*20)
                print('TF%ix%i vs. TF%ix%i'%(i_pt,i_rho,j_pt,j_rho))
                if((j_pt+1)*(j_rho+1) >= ft._nbins):
                    print('Number of parameters on alternative Bernstein configuration is exceeding number of bins. skipping')
                    continue
                if((j_pt+1)*(j_rho+1) == (i_pt+1)*(i_rho+1)):
                    print('Number of parameters is equal in both Bernstein configuration. skipping')
                    continue
                ft.orders_1 = (i_pt,i_rho)
                ft.orders_2 = (j_pt,j_rho)
                success = ft.process_files(dir_pattern,model_dir_pattern,skip_failed_fits = True)
                if(success>=0):
                    ft.plot_ftest()
                if(success<=0):
                    ft.plot_gofs()
                    
    ft.print_results()
    
if(__name__=='__main__'):

    

    # for algo in ["AD","KS","saturated"]:
    for algo in ["saturated","KS"]:
        process_TF_prep('/nfs/dust/cms/user/albrechs/JetMassCalibration/FTest_QCDTFScan_%s_Seed*'%algo,
                        'VJetsSelectionMCTFPt%iRho%iDataResTFPt0Rho0/qcdmodel/',
                        "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/QCDTfFTestPlots%s/"%algo,
                        algo)

        if algo == 'saturated':
        
            process_TF_prep('/nfs/dust/cms/user/albrechs/JetMassCalibration/FTest_DataTFScan_QCDOrder3x1_%s_Seed*'%algo,
                            'VJetsSelectionMCTFPt3Rho1DataResTFPt%iRho%i/',
                            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/Data2TfFTestPlots%s/"%algo,
                            algo,
                            n_p_other = 1 + (3+1)*(1+1))
        

            process_TF_prep('/nfs/dust/cms/user/albrechs/JetMassCalibration/FTest_DataTFScan_%s_Seed*'%algo,
                            'VJetsSelectionTFPt%iRho%i/',
                            "/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2_17/CMSSW_10_2_17/src/UHH2/JetMass/rhalph/DataTfFTestPlots%s/"%algo,
                            algo)
