#!/usr/bin/env python
from __future__ import print_function
import os,glob,sys
for i,p in enumerate(sys.path):
    if('sframebatch' in p.lower()):
        sys.path.pop(i)
sys.path.append(os.getenv("CMSSW_BASE")+'/src/UHH2/scripts/crab')
from readaMCatNloEntries import *

snippet='<In FileName="%s" Lumi="0.0"/>\n'
if(len(sys.argv)<2):
    exit(0)
inXML_path=os.path.abspath(sys.argv[1])
OutDir=inXML_path.replace(inXML_path.split('/')[-1],'PreSel/')
if(not os.path.exists(OutDir)):
    os.mkdir(OutDir)
inXML=open(inXML_path,'r+')
dir=''
for l in inXML:
    if('OUTdir' in l or 'OutputDirectory' in l):
        dir = l.split('"')[1]
        break
datasets=[]
for l in inXML:
    if('workdir' in l):
        dir+=l.split('"')[9]
    if('<InputData Version=' in l and '<!--' not in l):
        datasets.append(l.split('"')[1])
inXML.close()

for dataset in datasets:
    rule=dir+'/*'+dataset+'*'
    print rule
    files=glob.glob(rule)
    outXML=open(OutDir+dataset+'.xml','w')
    for file in files:
        outXML.write(snippet%file)
    outXML.close()
    readEntries(4,[OutDir+dataset+'.xml'],True)
