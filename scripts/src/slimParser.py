import numpy as np
import pandas as pd
import datatable as dt
import subprocess
import random
from string import Template
import re
from tempfile import NamedTemporaryFile
from multiprocessing import Pool
import glob
from tqdm import tqdm
import string

def runSlim(recipe,wStrength,sStrength,N,pposL,pposH,bgsMutRate,output,iteration,slimPath="/home/jmurga/.conda/envs/abc-mk/bin/slim"):
	"""
	Runs SLiM using subprocess.
	Args:
		command_line_args: list; a list of command line arguments.
	return: The entire SLiM output as a string.
	"""

	slimRecipe = Template(open(recipe, "r").read())

	mapping = {
		'weaklyStrength' : wStrength,
		'strongStrength' : sStrength,
		'N' : N,
		'pposL' : pposL,
		'pposH' : pposH,
		'bgsMutationRate' : bgsMutRate,
		
	}

	#print(slimRecipe.substitute(mapping))

	with NamedTemporaryFile("w") as slim_file:
		print(slimRecipe.substitute(mapping), file= slim_file,flush=True)

		slim = subprocess.run([slimPath, "-s", str(random.randint(1, 10**13)),slim_file.name],universal_newlines=True,stdout=subprocess.PIPE)

		slimResults = slim.stdout.split('\n')
		slimDaf = slimResults[slimResults.index('daf\tpi\tp0\tpw'):slimResults.index('di\td0\tdw')]
		slimDiv = slimResults[slimResults.index('di\td0\tdw'):slimResults.index('trueAlphaW\ttrueAlphaS\ttrueAlpha')]
		slimAlpha = slimResults[slimResults.index('trueAlphaW\ttrueAlphaS\ttrueAlpha'):-1]

		# Extract daf info from slim results
		slimDaf = [x.split('\t') for x in slimDaf]

		# Header daf
		h = slimDaf[0]
		# Data
		d = slimDaf[1:]
		# Pandas dataframe
		daf = pd.DataFrame(d,columns=h)
		daf.daf = np.around(daf.daf.astype(float),4)

		# Extract divergence info from slim results
		slimDiv = [x.split('\t') for x in slimDiv]
		slimAlpha = [x.split('\t') for x in slimAlpha]

		# Header div
		h = slimDiv[0]
		# Data
		d = slimDiv[1]
		# Pandas dataframe
		div = pd.DataFrame(np.reshape(np.array(d),(1,3)),columns=h)

		# Header alpha
		h = slimAlpha[0]
		# data
		d = slimAlpha[1]
		alpha = pd.DataFrame(np.reshape(np.array(d),(1,3)),columns=h)

		# Pandas dataframe
		daf.to_csv(output + "/daf/daf" + str(iteration) + ".tsv.gz",index=False,sep="\t",compression="gzip")
		div.to_csv(output + "/div/div" + str(iteration) + ".tsv.gz",index=False,sep="\t",compression="gzip")
		alpha.to_csv(output + "/div/alpha" + str(iteration) + ".tsv.gz",index=False,sep="\t",compression="gzip")

	# Parsing string output, we checked position on slim custom printed output and procces each variable taking into account correspondent positions. Excluding recipe execution info

def inputSlim(iterations,recipe,wStrength,sStrength,N,pposL,pposH,bgsMutRate,output,slimPath="/home/jmurga/.conda/envs/abc-mk/bin/slim"):

	totalSim = iterations[1]-iterations[0] +1
	irecipe = [recipe] * totalSim
	iwStrength = [wStrength] * totalSim
	isStrength = [sStrength] * totalSim
	iN = [N] * totalSim
	ipposL = [pposL] * totalSim
	ipposH = [pposH] * totalSim
	ibgsMutRate = [bgsMutRate] * totalSim
	ioutput = [output] * totalSim
	islimPath = [slimPath] * totalSim

	total = [i for i in range(iterations[0],iterations[1] + 1)]

	inp = list(zip(irecipe,iwStrength,isStrength,iN,ipposL,ipposH,ibgsMutRate,ioutput,total,islimPath))

	return(inp)

def parsePolDiv(path,N):
	"""
	slr, tupple array with daf and div data by element in list
	"""

	dafFiles = glob.glob(path + "/daf/*.tsv.gz")
	divFiles = glob.glob(path + "/div/div*.tsv.gz")
	alFiles  = glob.glob(path + "/div/al*.tsv.gz")

	iteration = len(dafFiles)

	sfs       = np.array(np.zeros(((N*2)-1,4)))
	divs      = np.array(np.zeros((1,3)))
	alphas    = np.array(np.zeros((iteration,3)))

	for f in tqdm(range(1,iteration )):
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
	divs      = dt.Frame({'di':divs[0],'d0':divs[1],'dw':divs[2]})
	alphas    = dt.Frame(alphas,names=['trueAlphaW','trueAlphaS','trueAlpha'])

	return(sfs.to_pandas(),divs.to_pandas(),alphas.to_pandas())

def parseSimulations(table,N):
	out = []; alphas = []

	for index,row in table.iterrows():
		try:
			sfs, div ,alphas = parsePolDiv(row.path,N)
			sfs.to_csv(row.path + "/sfs.tsv",header=True,index=False,sep="\t")
			div.to_csv(row.path + "/div.tsv",header=True,index=False,sep="\t")
			alphas.to_csv(row.path + "/alphas.tsv",header=True,index=False,sep="\t")
					
			try:
				cSfs = cumulativeSfs(sfs.to_numpy())
				asymp1 = amkt(cSfs[:,[0,4,2]],div.to_numpy().flatten(),0,1)
				asymp2 = amkt(cSfs[:,[0,1,2]],div.to_numpy().flatten(),0,1)
			except:
				cSfs = cumulativeSfs(reduceSfs(sfs.to_numpy(),100))
				asymp1 = amkt(cSfs[:,[0,4,2]],div.to_numpy().flatten(),0,1)
				asymp2 = amkt(cSfs[:,[0,1,2]],div.to_numpy().flatten(),0,1)

			tmp = pd.DataFrame(np.mean(alphas,axis=0)).T
			tmp['asymp'] = asymp[0]['alpha']
			tmp['analyticalEstimation']= row.estimation
			tmp['alpha_nopos']= asymp1[1]
			tmp['alpha']= asymp2[1]
			al.append(tmp)
		except:
			continue
	df = pd.concat(out)
	df = df.set_index(table.path.apply(lambda row: row.split('/')[-1]))
	return(df)

def randomString():
    return ''.join(random.choice(string.ascii_letters) for m in range(0,8))

def priorsJulia(table,nSimulations,model,script="/home/jmurga/mkt/202004/scripts/src/sim.jl",regfile="job.txt",threads=4):
	
	for index,row in table.iterrows():
		print(index, row.path)
		nScript    = [script] * threads
		nNames     = [randomString() for i in range(0,threads)]
		nModel     = [model] * threads

		nResults   = [str(nSimulations)] * threads
		nPaths     = [row.path.split('/')[-1]] * threads
		
		tmp        = [list(x) for x in zip(nScript,nModel,nPaths,nResults,nNames)]

		# Create reg file
		jobFile=re.sub('simulations','summStat',row.path) + '/' + regfile 
		output=re.sub('simulations','summStat',row.path)
		

		for i in nNames:
			p =  output + '/' + output.split('/')[-1] + '_' + i + '1.tsv' 
			d =  output + '/' + 'sfsNoDemog.tsv'
			
			f = open(jobFile, "a")
			f.write('/home/jmurga/ABCreg/src/reg -p ' + p + ' -d ' + d + ' -P 3 -S 104 -b ' + output + '/' + i + ' -T -t 0.001' + '\n')
			f.close()


		p = Pool(processes=threads)
		output = p.map(runJulia,tmp)
		p.terminate()

def runJulia(simulatedList,juliaPath="/home/jmurga/julia-1.5.0/bin/julia"):
	subprocess.run([juliaPath] + ["--sysimage"] + ["/home/jmurga/mktest.so"] +simulatedList)

