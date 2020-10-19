import sys
sys.path.insert(0, '/home/jmurga/mkt/202004/scripts/src/pyAmkt.py')
from pyAmkt import *
import numpy as np
import pandas as pd
import datatable as dt
import subprocess
import random
from string import Template
import re
import os
from tempfile import NamedTemporaryFile
from multiprocessing import Pool
import glob
from tqdm import tqdm
import string


def runSlim(recipe,simTable,pSize,codingLength,strongStrength,weaklyStrength,bins,replicas,threads,slimPath="/home/jmurga/.conda/envs/abcmk/bin/slim",parallelPath="/home/jmurga/.conda/envs/abcmk/bin/parallel"):

    bgs = []
    os.makedirs(simTable.path.values[0],exist_ok=True)    
    os.makedirs(simTable.path.values[0] + '/daf',exist_ok=True)    
    os.makedirs(simTable.path.values[0] + '/div',exist_ok=True)    
    os.makedirs(simTable.path.values[0] + '/vcf',exist_ok=True)    

    path =  "\\" + "'" + simTable.path.values[0] + "\\" + "'" 
    tmp = parallelPath + " -j" + str(threads) + " \"" + slimPath + " -d popSize=" +str(pSize) + " -d codingLength=" + str(codingLength) + " -d bins=" + str(bins) + " -d weaklyStrength=" + str(weaklyStrength)  + " -d strongStrength=" + str(500)  + " -d bgsMutationRate=" + str(simTable.bgsThetaF.values[0])  + " -d pposL=" + str(simTable.pposL.values[0]) + " -d pposH=" + str(simTable.pposH.values[0]) + " -d p=" + path + " -d nF={} " + recipe + "\" ::: {"+ str(replicas[0]) + ".."+ str(replicas[1]) + "}"
    print(tmp)

    process = subprocess.run(tmp, shell=True,check=True,executable='/bin/bash')

def parsePolDiv(path,N,sample=None):
    """
    slr, tupple array with daf and div data by element in list
    """

    dafFiles = np.sort(glob.glob(path + "/daf/*.tsv.gz"))
    divFiles = np.sort(glob.glob(path + "/div/div*.tsv.gz"))
    alFiles  = np.sort(glob.glob(path + "/div/al*.tsv.gz"))

    if sample is not None:
        tmp = np.arange(0,len(dafFiles))
        idx   = np.sort(np.random.choice(tmp,sample,replace=False))
        dafFiles = dafFiles[idx]
        divFiles = divFiles[idx]
        alFiles  = alFiles[idx]

        iteration= len(dafFiles)
    else:
        iteration = len(dafFiles)
    

    sfs       = np.array(np.zeros((N,4)))
    divs      = np.array(np.zeros((1,4)))
    alphas    = np.array(np.zeros((iteration,3)))

    for f in tqdm(range(0,iteration)):
        daf             = dt.fread(dafFiles[f], sep='\t',columns=dt.float64)
        daf[:,'pi_nopos'] = daf[:, dt.f.pi - dt.f.pw]
        daf = daf.to_numpy()
        tmp = daf[:,1:]
        sfs = sfs + tmp

        d = dt.fread(divFiles[f], sep='\t',columns=dt.float64).to_numpy()
        divs = divs + d

        al = dt.fread(alFiles[f], sep='\t',columns=dt.float64).to_numpy()
        alphas[f,:] = al
        # daf             = pd.read_csv(dafFiles[f], header=0, sep='\t')
        # daf['pi_nopos'] = daf.pi - daf.pw
        # sfs             = sfs + daf.iloc[:,1:]

        # d = pd.read_csv(divFiles[f], header=0, sep='\t')
        # divs = divs + d

        # al = pd.read_csv(alFiles[f], header=0, sep='\t')
        # alphas.iloc[f] = al.to_numpy()

    sfs = np.hstack([daf[:,0].reshape(daf.shape[0],1),sfs])
    sfs = dt.Frame(sfs,names=['f','pi','p0','pw','pi_nopos'])
    
    divs = np.sum(divs,axis=0)
    divs      = dt.Frame({'di':divs[0],'d0':divs[1],'dw':divs[2],'ds':divs[3]})
    alphas    = dt.Frame(alphas,names=['trueAlphaW','trueAlphaS','trueAlpha'])

    if sample is not None:
        sfs.to_pandas().to_csv(path + "/sfs" + str(sample) + ".tsv",header=True,index=False,sep="\t")
        divs.to_pandas().to_csv(path + "/div" + str(sample) + ".tsv",header=True,index=False,sep="\t")
        alphas.to_pandas().to_csv(path + "/alphas"+ str(sample) + ".tsv",header=True,index=False,sep="\t")
    else:
        sfs.to_pandas().to_csv(path + "/sfs.tsv",header=True,index=False,sep="\t")
        divs.to_pandas().to_csv(path + "/div.tsv",header=True,index=False,sep="\t")
        alphas.to_pandas().to_csv(path + "/alphas.tsv",header=True,index=False,sep="\t")
    
def saveSimulatedAlphas(table,bins,reduced=None,sample=None):
    out = [];al = [];
    for index,row in table.iterrows():
        print(row.path)

        if sample is not None:
            sfs = pd.read_csv(row.path + "/sfs"+str(sample)+".tsv",header=0,sep="\t")
            div = pd.read_csv(row.path + "/div"+str(sample)+".tsv",header=0,sep="\t")
            alphas = pd.read_csv(row.path + "/alphas"+str(sample)+".tsv",header=0,sep="\t")
            alpha = (div.ds+div.dw)/div.di
        
        else:
            sfs = pd.read_csv(row.path + "/sfs.tsv",header=0,sep="\t")
            div = pd.read_csv(row.path + "/div.tsv",header=0,sep="\t")
            alphas = pd.read_csv(row.path + "/alphas.tsv",header=0,sep="\t")
            alpha = (div.ds+div.dw)/div.di
        
        if(reduced is not None):
            cSfs = cumulativeSfs(sfs.to_numpy())
            asymp1 = amkt(cSfs[:,[0,4,2]],div.to_numpy().flatten(),0,1)
            asymp2 = amkt(cSfs[:,[0,1,2]],div.to_numpy().flatten(),0,1)
            f = np.arange(1,sfs.shape[0]+1)
        else:
            cSfs = cumulativeSfs(reduceSfs(sfs.to_numpy(),bins))
            asymp1 = amkt(cSfs[:,[0,4,2]],div.to_numpy().flatten(),0,1)
            asymp2 = amkt(cSfs[:,[0,1,2]],div.to_numpy().flatten(),0,1)
            f = np.arange(1,bins+1)
            
        tmp = pd.DataFrame({'trueAlpha':alpha,'asymp_nopos':asymp1[0]['alpha'],'asymp':asymp2[0]['alpha'],'analyticalEstimation':row.estimation,'path':row.path.split('/')[-1]})
        
        tmpAlpha = pd.DataFrame({'f':f,'alphas':asymp1[1],'sfs':np.repeat(['nopos'],f.shape[0]),'path':row.path.split('/')[-1]})
        al.append(tmpAlpha)
        tmpAlpha = pd.concat([tmpAlpha,pd.DataFrame({'f':f,'alphas':asymp2[1],'sfs':np.repeat(['pos'],f.shape[0]),'path':row.path.split('/')[-1]})])
        al.append(tmpAlpha)
        out.append(tmp)
    return(pd.concat(out),pd.concat(al))

def randomString():
    return ''.join(random.choice(string.ascii_letters) for m in range(0,8))

def priorsJulia(table,nSimulations,pSize,nSize,model,bins,script="/home/jmurga/mkt/202004/scripts/src/sim.jl",regfile="rJob.sh",threads=4,replicas=[1,4],precomipledImg="/home/jmurga/mkt/202004/scripts/src/mktest.so",parallelPath="/home/jmurga/.conda/envs/abcmk/bin/parallel",abcreg="/home/jmurga/ABCreg/src/reg"):
    
    for index,row in table.iterrows():
        # Create reg file

        output=re.sub('simulations','summStat',row.path)
        os.makedirs(output,exist_ok=True)

        if precomipledImg is not None:
            jl = "julia -J " + precomipledImg
        else:
            jl = "julia "
        
        tmp = parallelPath + " -j" + str(threads) + " \"" + jl + " " + script + " " + model + " " + row.path.split('/')[-1] + " " + str(pSize) + " " + str(nSize) + " " + str(bins) + " " + str(nSimulations) + " {}" + "\" ::: {"+ str(replicas[0]) + ".."+ str(replicas[1]) + "}"
        
        print(tmp)
        process = subprocess.run(tmp, shell=True,check=True,executable='/bin/bash')
        
        #Merge files
        # tmpMerge = "cat " + output + "/" + row.path.split('/')[-1] + "_* > " + output + '/' "prior_" + str(bins) + ".tsv"
        # print(tmpMerge)
        # process = subprocess.run(tmpMerge, shell=True,check=True,executable='/bin/bash')
        
        # #Delete files
        # tmpDelete = "rm " + output +  "/" + row.path.split('/')[-1] + "_*"
        # # print(tmpDelete)
        # process = subprocess.run(tmpDelete, shell=True,check=True,executable='/bin/bash')

        p =  output + '/' + row.path.split('/')[-1] + '.tsv' 
        d =  output + '/' + 'sfs' + model + '.tsv'
        
        f = open(output + "/" + regfile, "a")
        f.write(abcreg + ' -p ' + p + ' -d ' + d + ' -P 3 -S 100 -b ' + output + '/computation' + ' -T -t 0.001' + '\n')
        f.close()


        # for i in nNames:
        #     p =  output + '/' + output.split('/')[-1] + '_' + i + '1.tsv' 
        #     d =  output + '/' + 'sfsNoDemog.tsv'
            
        #     f = open(jobFile, "a")
        #     f.write('/home/jmurga/ABCreg/src/reg -p ' + p + ' -d ' + d + ' -P 3 -S 100 -b ' + output + '/' + i + ' -T -t 0.001' + '\n')
        #     f.close()


        # p = Pool(processes=threads)
        # output = p.map(runJulia,tmp)
        # p.terminate()

def abcOutputs(model,simulations):
    
    alphas = dict.fromkeys(simulations)
    density = dict.fromkeys(simulations)
    plots = dict.fromkeys(simulations)
    
    for n in simulations:

        print(n)
        sim= PATH + '/rawData/simulations/' + model + '/' + n
        ss= PATH + '/rawData/summStat/' + model + '/' + n
        # alphas[[n]] = fread(paste0(sim,"/alphas.tsv")) %>% summarize_all(mean)
        tmp = pd.read_csv(sim + "/alphas.tsv",sep='\t') 
        tmp[['analysis']] = n
        tmp[['method']] = "simulation"
        tmp.columns = ['alphaW','alphaS','alpha','analysis',"method"]

        abcResults = glob.glob(ss + "/*.tangent*")
        l = list()
        d = list()

        for i in abcResults:
            df = pd.read_csv(i,sep='\t',names = ['alphaW','alphaS','alpha'])
            df[['analysis']] = n
            d.append(df)
            l.append(np.mean(df).to_numpy())
        
        density[n] = pd.concat(d).reset_index(drop=True)
        alphas[n] = pd.DataFrame(np.vstack(l),columns=['alphaW','alphaS','alpha'])

        alphas[n][['analysis']] = n
        alphas[n][['method']] = "abc"
        alphas[n] = pd.concat([alphas[n],tmp])
        
    alphaToPlot = pd.concat(alphas.values())
    densityToPlot = pd.concat(density.values())

    return(alphaToPlot,densityToPlot)