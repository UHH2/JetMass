import rhalphalib as rl
import numpy as np
from ROOT import TFile, TH1F
def uhh_producer():
    histLocation='/afs/desy.de/user/a/albrechs/xxl/af-cms/dazsle/JetMass/Histograms/top'
    Variable='Mass'
    HistDir='JetMass_pt200'
    ModelName='TopMass'
    Samples=['SingleTop','TTbar','WJets','other']

    ptBins=set()
    etaBins=set()
    file=TFile('%s/%s.root'%(histLocation,Samples[0]),'READ')
    for TKey in file.GetDirectory(HistDir).GetListOfKeys():
        key=TKey.GetName()
        if Variable in key:
            if'central' not in key:
                ptBins.add(int(key.split('_')[1]))
                etaBins.add(int(key.split('_')[2]))
    file.Close()

    model = rl.Model(ModelName)

    gridNuisances=[]
    for i in ptBins:
        for j in etaBins:
            gridNuisances.append([rl.NuisanceParameter('massScale_pt%i_eta%i'%(i,j),'shape'),i,j])

    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')

    ch=rl.Channel('ch1')
    model.addChannel(ch)
    for sName in Samples:
        sampleType=rl.Sample.SIGNAL if sName=='TTbar' else rl.Sample.BACKGROUND
        sFile=TFile('%s/%s.root'%(histLocation,sName),'READ')
        hist=sFile.Get(HistDir+'/'+Variable+'_central')
        sample=rl.TemplateSample(ch.name+'_'+sName,sampleType,hist)
        for (gridNuisance,pt,eta) in gridNuisances:
            histUp=sFile.Get('%s/%s_%i_%i_up'%(HistDir,Variable,pt,eta))
            histDown=sFile.Get('%s/%s_%i_%i_down'%(HistDir,Variable,pt,eta))
            sample.setParamEffect(gridNuisance,histUp,histDown)
        sample.setParamEffect(lumi, 1.027)
        ch.addSample(sample)

    dataFile=TFile('%s/%s.root'%(histLocation,'Data'))
    print(dataFile)
    ch.setObservation(dataFile.Get(HistDir+'/'+'Mass_central'))

    import pickle
    with open("model.pkl", "wb") as fout:
        pickle.dump(model, fout)

    import sys
    print("ROOT used? ", 'ROOT' in sys.modules)
    model.renderCombine(ModelName)
    print("ROOT used? ", 'ROOT' in sys.modules)


if __name__ == '__main__':
    uhh_producer()
