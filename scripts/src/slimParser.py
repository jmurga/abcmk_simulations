import sys
sys.path.insert(0, '/home/jmurga/mkt/202004/scripts/src/pyAmkt.py')
from pyAmkt import *
import numpy as np
import pandas as pd
# import datatable as dt
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

def saveSimulatedAlphas(table,l=0,h=1,bins=None,sample=None,cumulative=True):
	out = [];al = [];
	for index,row in table.iterrows():
		# print(row.path)
		r = np.random.choice(np.arange(1,100),1)[0]

		if sample is not None:

			sfs = pd.read_csv(row.path + "/sfs"+str(sample)+".tsv",header=0,sep="\t")
			div = pd.read_csv(row.path + "/div"+str(sample)+".tsv",header=0,sep="\t")
			# alphas = pd.read_csv(row.path + "/alphas"+str(sample)+".tsv",header=0,sep="\t")
			alpha = (div.ds+div.dw)/div.di
		
		else:
			
			sfs = pd.read_csv(row.path + "/sfs"+ str(r)+".tsv",header=0,sep="\t")
			div = pd.read_csv(row.path + "/div"+ str(r)+".tsv",header=0,sep="\t")
			# alphas = pd.read_csv(row.path + "/alphas.tsv",header=0,sep="\t")
			alpha = (div.ds+div.dw)/div.di
		
		if(bins is not None):
			if cumulative:
				cSfs = reduceSfs(cumulativeSfs(sfs.to_numpy()),bins)
			else:
				cSfs = reduceSfs(sfs.to_numpy(),bins)

			cSfs = cSfs[(cSfs[:,0] >= l) & (cSfs[:,0] <=h)]
			asymp1 = amkt(cSfs[:,[0,4,2]],div.to_numpy().flatten(),0,1)
			asymp2 = amkt(cSfs[:,[0,1,2]],div.to_numpy().flatten(),0,1)

			f = np.arange(1,asymp1[1].shape[0]+1)/bins

		else:
			if cumulative:
				cSfs = cumulativeSfs(sfs.to_numpy())
			else:
				cSfs = sfs.to_numpy()
				
			# cSfs = cSfs[(cSfs[:,0] >= l) & (cSfs[:,0] <=h)]

			asymp1 = amkt(cSfs[:,[0,4,2]],div.to_numpy().flatten(),l,h)

			asymp2 = amkt(cSfs[:,[0,1,2]],div.to_numpy().flatten(),l,h)	
			f = np.arange(1,asymp1[1].shape[0]+1)/cSfs.shape[0]
		
		tmp = pd.DataFrame({'trueAlpha':alpha,'asymp_nopos':asymp1[0]['alpha'],'asymp':asymp2[0]['alpha'],'analyticalEstimation':row.estimation,'path':row.path.split('/')[-1],'analysis':row.analysis})
		
		tmpAlpha = pd.DataFrame({'f':f,'alphas':asymp1[1],'B':np.repeat(row.B,f.shape[0]),'alphaW':np.repeat(row.alphaW,f.shape[0]),'sfs':np.repeat(['Neutral + deleterious'],f.shape[0]),'path':row.path.split('/')[-1],'analysis':row.analysis})

		tmpAlpha = pd.concat([tmpAlpha, pd.DataFrame({'f': f, 'alphas': asymp2[1], 'sfs':np.repeat(['All alleles'], f.shape[0]), 'B':np.repeat(
			row.B, f.shape[0]), 'alphaW':np.repeat(row.alphaW, f.shape[0]), 'path':row.path.split('/')[-1], 'analysis':row.analysis})])
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

def parseBootstrapPolDiv(path,N,sample=0.01,replicas=100,nthreads=8,dofe=False,bins=None):
	"""
	slr, tupple array with daf and div data by element in list
	"""
	dafFiles = np.sort(glob.glob(path + "/daf/*.tsv.gz"))
	divFiles = np.sort(glob.glob(path + "/div/div*.tsv.gz"))
	pool = Pool(processes = nthreads)
	lSfs,ldiv = zip(*pool.starmap(openFiles,list(zip(dafFiles,divFiles))))
	pool.terminate()

	# for i in range(1,divFiles.shape[0]):
	# 	print(i)
	# 	a,b = openFiles(dafFiles[i],divFiles[i])
	idxAll = np.arange(0,len(dafFiles)-1)
	bn = (N*2) - 1
	rsample = int(sample * len(dafFiles))


	sfs = np.sum(lSfs,axis=0)[:,1:]
	f = np.reshape(np.arange(1,bn+1)/bn,(bn,1))
	sfs = np.hstack((f,sfs))
	# sCumu = cumulativeSfs(sfs)
	d = np.sum(ldiv,axis=0)

	fName = path + path.split("/")[-1] 
	sfs = pd.DataFrame(np.round(sfs,5),columns=['f','pi','p0','pw'])
	sfs['pi_nopos'] = sfs['pi'] - sfs['pw']
	d = pd.DataFrame(d).T
	d.columns=['di','d0','dw','ds']

	if bins is not None:
		sfs.to_csv(fName + "/sfs.tsv", header=True, index=False, sep="\t")
		d.to_csv(fName + "/div.tsv", header=True, index=False, sep="\t")
	else:
		sfs.to_csv(path + "/sfs.tsv", header=True, index=False, sep="\t")
		d.to_csv(path + "/div.tsv", header=True, index=False, sep="\t")

	for r in tqdm(range(1, replicas+1)):
		idx   = np.sort(np.random.choice(idxAll,rsample,replace=True))
		tmpDaf = [lSfs[i] for i in idx]
		tmpDiv = [ldiv[i] for i in idx]
		m = np.array([rsample*2000])

		sfs = np.sum(tmpDaf,axis=0)[:,1:]
		f = np.reshape(np.arange(1,bn+1)/bn,(bn,1))
		sfs = np.hstack((f,sfs))
		# sCumu = cumulativeSfs(sfs)
		d = np.sum(tmpDiv,axis=0)
		
		if dofe is True:
			os.makedirs(path + path.split("/")[-1], exist_ok=True) 
			if bins is not None:
				header = pd.DataFrame([1,1,bins*2]).T
				fName = path + "/" + path.split("/")[-1] + "/" + path.split("/")[-1] + "_polydfe_"+ str(r) + "_"  + str(bins) + ".tsv"
				sfs = reduceSfs(sfs,bins)
			else:
				header = pd.DataFrame([1,1,N*2]).T
				fName = path + "/" + path.split("/")[-1] + "/" + path.split("/")[-1] + "_polydfe_"+ str(r) + ".tsv"

			neutral = pd.DataFrame(np.hstack([sfs[:,2],m*0.25,d[1],m*0.25])).T
			selected = pd.DataFrame(np.hstack([sfs[:,1],m*0.75,d[0],m*0.75])).T
			df = pd.concat([neutral,selected])
			header.to_csv(fName,sep=' ',header=None,index=False)
			df.to_csv(fName,sep=' ',header=None,index=False,mode='a')			
		else:
			fName = path + path.split("/")[-1] 
			sfs = pd.DataFrame(np.round(sfs,5),columns=['f','pi','p0','pw'])
			sfs['pi_nopos'] = sfs['pi'] - sfs['pw']
			d = pd.DataFrame(d).T
			d.columns=['di','d0','dw','ds']
			if bins is not None:
				sfs.to_csv(fName + "/sfs_" + str(r) + ".tsv", header=True, index=False, sep="\t")
				d.to_csv(fName + "/div_" + str(r) + ".tsv", header=True, index=False, sep="\t")
			else:
				sfs.to_csv(path + "/sfs" + str(r) + ".tsv", header=True, index=False, sep="\t")
				d.to_csv(path + "/div" + str(r) + ".tsv", header=True, index=False, sep="\t")

def bootstrap(ls,ld,lIndx,bn,smpl,path,r,output,dofe,bins):
	idx   = np.sort(np.random.choice(lIndx,smpl,replace=True))
	tmpDaf = [ls[i] for i in idx]
	tmpDiv = [ld[i] for i in idx]
	m = np.array([smpl*2000])

	sfs = np.sum(tmpDaf,axis=0)[:,1:]
	f = np.reshape(np.arange(1,bn+1)/bn,(bn,1))
	sfs = np.hstack((f,sfs))
	# sCumu = cumulativeSfs(sfs)
	d = np.sum(tmpDiv,axis=0)
	
	if dofe is True:
		os.makedirs(output + path.split("/")[-1], exist_ok=True) 
		if bins is not None:
			header = pd.DataFrame([1,1,bins*2]).T
			fName = output + "/" + path.split("/")[-1] + "/" + path.split("/")[-1] + "_polydfe_"+ str(r) + "_"  + str(bins) + ".tsv"
			sfs = reduceSfs(sfs,bins)
		else:
			header = pd.DataFrame([1,1,N*2]).T
			fName = output + "/" + path.split("/")[-1] + "/" + path.split("/")[-1] + "_polydfe_"+ str(r) + ".tsv"

		neutral = pd.DataFrame(np.hstack([sfs[:,2],m*0.25,d[1],m*0.25])).T
		selected = pd.DataFrame(np.hstack([sfs[:,1],m*0.75,d[0],m*0.75])).T
		df = pd.concat([neutral,selected])
		header.to_csv(fName,sep=' ',header=None,index=False)
		df.to_csv(fName,sep=' ',header=None,index=False,mode='a')	
	else:
		fName = output + path.split("/")[-1] 
		sfs = pd.DataFrame(np.round(sfs,5),columns=['f','pi','p0','pw'])
		sfs['pi_nopos'] = sfs['pi'] - sfs['pw']
		d = pd.DataFrame(d).T
		d.columns=['di','d0','dw','ds']
	
		sfs.to_csv(path + "/sfs" + str(r) + ".tsv", header=True, index=False, sep="\t")
		d.to_csv(path + "/div" + str(r) + ".tsv", header=True, index=False, sep="\t")

def compareAlphas(simulations,nn,PATH="/home/jmurga/mkt/201902/rawData/simulations/"):
    
    df = pd.DataFrame(np.zeros((simulations.shape[0],6)),columns = ['analysis','amk','amk cutoff','imputed mk','imputed mk high frequencies','true alpha'])
    polData = []
    tmp = []

    for s in range(0,simulations.shape[0]):
        print(simulations[s])
        sfs = pd.read_csv(PATH + simulations[s] + '/sfs.tsv',sep='\t').to_numpy()
        div = pd.read_csv(PATH + simulations[s] + '/div.tsv',sep='\t')
        d = div.to_numpy()[0]

        sCumu = cumulativeSfs(sfs)


        amkFull = amkt(sCumu,d,0,1)
        amk = round(amkFull[0]['alpha'],3)
        amkCutoff = round(amkt(sCumu,d,0,0.9)[0]['alpha'],3)
        imk, data1 = imputedMKT(sfs,d,0.15)
        imkHigh, data2 = imputedMKT(sfs,d,0.15,0.9)

        if ('dw' in div.columns):
            trueAlpha = round((d[2]+d[3])/d[0],3)
        else:
            trueAlpha = round(d[2]/d[0],3)

        df.iloc[s,:] = np.array([simulations[s],amk,amkCutoff,imk,imkHigh,trueAlpha])

        polData.append(np.hstack([simulations[s],'[' +str(0.15)+'-'+str(1) + ']',data1]))
        polData.append(np.hstack([simulations[s],'[' +str(0.15)+'-'+str(0.9) + ']',data2]))

        tmp.append(pd.DataFrame({'f':np.arange(1,nn+1),'alpha':amkFull[1],'analysis':[simulations[s]]*nn}))

    dfAlpha = pd.concat(tmp)
    polData = pd.DataFrame(polData,columns = ['analysis','cutoff','pi', 'piNeutral', 'p0', 'piLow/p0Low', 'piInter/p0Inter', 'piHigh/p0High', 'di', 'd0'])
    return df,dfAlpha,polData

def sfsToDofe(path,bins,output="/home/jmurga/mkt/202004/rawData/dofe/grapes/"):


	sfs= pd.read_csv(path + "/sfs.tsv",delimiter='\t').to_numpy()
	div= pd.read_csv(path + "/div.tsv",delimiter='\t').to_numpy().flatten().astype(int)

	if ((bins*2)-1) != sfs.shape[0]:
		sfs = reduceSfs(sfs,bins)
		header = pd.DataFrame(["#unfolded"])
		fName = output + path.split("/")[-1] + "_grapes_"+ str(bins) + ".tsv"
	else:
		header = pd.DataFrame(["#unfolded"])
		fName = output + path.split("/")[-1] + "_grapes.tsv"

	sfs = sfs.astype(int)
	data = pd.DataFrame(np.hstack([path.split("/")[-1],bins*2,int(10**8*0.75),sfs[:,1],int(10**8*0.25),sfs[:,2],int(10**8*0.75),div[0],int(10**8*0.25),div[1]])).T

	f = open(fName,"w")
	f.write(" \n")
	f.close()
	header.to_csv(fName,sep='\t',header=None,index=False)
	data.to_csv(fName,sep='\t',header=None,index=False,mode='a')

def openFiles(sfsFile,divFile):

	sfs = pd.read_csv(sfsFile,sep='\t').to_numpy()
	div = pd.read_csv(divFile,sep='\t').to_numpy().flatten()

	return(sfs,div)

def dofeToSfs(file,output):

	# file = "primates_fruitflies.dofe"

	f = open(file)
	content = f.readlines()

	content = np.array(content[2:])
	content = content[content!='\n']
	for s in tqdm(content):
		tmp = s.split('\t')
		header = tmp[0].lower()
		samples = int(tmp[1])
		data = np.array(tmp[2:])
		pn = data[0:samples][1:].astype(float)
		ps = data[samples:(samples*2)][1:].astype(float)

		d = data[samples*2:][[1,3]].astype(float)

		sfs = pd.DataFrame({'f':np.round(np.arange(1,samples)/samples,3),'pn':pn,'ps':ps})
		divergence = pd.DataFrame({'dn':d[0],'ds':d[1]},index=[0])
		os.makedirs(output + "/" + header ,exist_ok=True)    
		sfs.to_csv(output + "/" + header + "/sfs_" + header + ".tsv",header=True,index=False,sep='\t')
		divergence.to_csv(output + "/" + header + "/div_" + header  + ".tsv",header=True,index=False,sep='\t')

		f = open(output + "/" + header + "/" + header  + ".dofe","w")
		f.write('#unfolded\n'+s)
		f.close()