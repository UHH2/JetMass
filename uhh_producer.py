import rhalphalib as rl
import numpy as np
import sys
from ROOT import TFile, TH1F
def uhh_producer(channelConfigs=None,ModelName='UHH_Model',gridHistFileName='/afs/desy.de/user/a/albrechs/xxl/af-cms/UHH2/10_2/CMSSW_10_2_10/src/UHH2/JetMass/Histograms/grid.root'):

    print('reading grid for nuisance parameter:')
    gridHistFile=TFile(gridHistFileName,'READ')
    gridHist=gridHistFile.Get('grid')
    gridAxises=dict(item.strip().split("=") for item in gridHist.GetTitle().split(","))
    xBins=range(gridHist.GetNbinsX())
    yBins=range(gridHist.GetNbinsY())

    categoriesHist=gridHistFile.Get('categories')
    ParticleCategories=[]
    for i in range(1,categoriesHist.GetNbinsX()+1):
        ParticleCategories.append(categoriesHist.GetXaxis().GetBinLabel(i))

    gridHistFile.Close()
    print('used variation categories:',ParticleCategories)
    print('X: %s , %i bins'%(gridAxises['x'],len(xBins)))
    print('Y: %s , %i bins'%(gridAxises['y'],len(yBins)))
    if(channelConfigs == None):
        print('must specify channel Configurations!')
        exit(0)

    model = rl.Model(ModelName)
    gridNuisances=[]
    print('adding nuisance parameters:')
    for category in ParticleCategories:
        for i in xBins:
            for j in yBins:
                print('massScale_%s%i_%s%i_%s'%(gridAxises['x'],i,gridAxises['y'],j,category),'shape')
                gridNuisances.append([rl.NuisanceParameter('massScale_%s%i_%s%i_%s'%(gridAxises['x'],i,gridAxises['y'],j,category),'shape'),i,j,category])


    print('CMS_lumi', 'lnN')
    lumi = rl.NuisanceParameter('CMS_lumi', 'lnN')

    for channelName,config in channelConfigs.items():
        histLocation=config['histLocation']
        Variable=config['variable']
        histDir=config['histDir']
        Samples=config['samples']

        if ('regions' not in config) or (len(config['regions'])==0):
            config['regions']=['']

        for region in config['regions']:
            regionName=channelName+region
            ch=rl.Channel(regionName)
            model.addChannel(ch)
            for sName in Samples:
                sampleType=rl.Sample.SIGNAL if sName==config['signal'] else rl.Sample.BACKGROUND
                sFile=TFile('%s/%s.root'%(histLocation,sName),'READ')
                hist=sFile.Get(histDir+'/'+Variable+'_central')
                sample=rl.TemplateSample(ch.name+'_'+sName,sampleType,hist)
                for (gridNuisance,x,y,category) in gridNuisances:
                    histUp=sFile.Get('%s/%s_%i_%i_%s_up'%(histDir,Variable,x,y,category))
                    histDown=sFile.Get('%s/%s_%i_%i_%s_down'%(histDir,Variable,x,y,category))
                    sample.setParamEffect(gridNuisance,histUp,histDown)
                    sample.setParamEffect(lumi, 1.027)
                ch.addSample(sample)

            dataFile=TFile('%s/%s.root'%(histLocation,config['obs']))

            if 'varyPseudoLike' in config:
                histPath=config['varyPseudoLike']
                if '/' not in histPath:
                    histPath=histDir+'/'+histPath
                dataHist=dataFile.Get(histPath)
                dataHist.SetName('Mass_central')
            else:
                dataHist=dataFile.Get(histDir+'/'+'Mass_central')

            ch.setObservation(dataHist)

    # import pickle
    # with open("model.pkl", "wb") as fout:
    #     pickle.dump(model, fout)

    print("ROOT used? ", 'ROOT' in sys.modules)
    model.renderCombine(ModelName)
    print("ROOT used? ", 'ROOT' in sys.modules)


if __name__ == '__main__':
    import json
    configs=json.load(open(sys.argv[1]))
    uhh_producer(configs)
    from runFit import runFits
    runFits([configs['ModelName']],configs['pathCMSSW'])
    # for modelName in ['UHH_Model_0','UHH_Model_1','UHH_Model_2']:
    #     uhh_producer(channels,ModelName=modelName)
