import sys
PATH = "/home/jmurga/mkt/202004" 
sys.path.insert(0, PATH + '/scripts/src/')  
from py_amkt import *
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

def save_simulated_alphas(table):
	out = [];al = [];
	for index,row in table.iterrows():
		
		sfs = pd.read_csv(row.path + "/sfs.tsv",header=0,sep="\t")
		div = pd.read_csv(row.path + "/div.tsv",header=0,sep="\t")
		# alphas = pd.read_csv(row.path + "/alphas.tsv",header=0,sep="\t")
		alpha = (div.ds+div.dw)/div.di
		c_sfs = cumulative_sfs(sfs.to_numpy())

		asymp1 = amkt(c_sfs[:,[0,4,2]],div.to_numpy().flatten(),0,1)

		asymp2 = amkt(c_sfs[:,[0,1,2]],div.to_numpy().flatten(),0,1)	
		f = np.arange(1,asymp1[1].shape[0]+1)/c_sfs.shape[0]
	
		tmp = pd.DataFrame({'trueAlpha':alpha,'asymp_nopos':asymp1[0]['alpha'],'asymp':asymp2[0]['alpha'],'analyticalEstimation':row.estimation,'path':row.path.split('/')[-1]})
		
		tmpAlpha = pd.DataFrame({'f':f,'alphas':asymp1[1],'B':np.repeat(row.B,f.shape[0]),'alpha_weak':np.repeat(row.alpha_weak,f.shape[0]),'sfs':np.repeat(['Neutral + deleterious'],f.shape[0]),'path':row.path.split('/')[-1]})

		tmpAlpha = pd.concat([tmpAlpha, pd.DataFrame({'f': f, 'alphas': asymp2[1], 'sfs':np.repeat(['All alleles'], f.shape[0]), 'B':np.repeat(
			row.B, f.shape[0]), 'alpha_weak':np.repeat(row.alpha_weak, f.shape[0]), 'path':row.path.split('/')[-1]})])
		al.append(tmpAlpha)
		out.append(tmp)
	return(pd.concat(out),pd.concat(al))

def parse_and_bootstrap(path,N,replicas=100,nthreads=8):
	"""
	slr, tupple array with daf and div data by element in list
	"""
	daf_files = np.sort(glob.glob(path + "/daf/*.tsv.gz"))
	div_files = np.sort(glob.glob(path + "/fix/div*.tsv.gz"))
	
	# pool = Pool(processes = nthreads)
	# l_sfs,l_div = zip(*pool.starmap(open_files,list(zip(daf_files,div_files))))
	# pool.terminate()

	with Pool(processes=nthreads) as pool: 
		l_sfs,l_div = zip(*pool.starmap(open_files,zip(daf_files,div_files)))

	idx_all = np.arange(0,len(daf_files)-1)
	bn      = (N*2)

	sfs     = np.sum(l_sfs,axis=0)[:,1:]
	f       = np.reshape(np.arange(1,bn)/bn,(bn-1,1))
	sfs     = np.hstack((f,sfs))
	d       = np.sum(l_div,axis=0)
	d       = pd.DataFrame(d).T
	
	try:
		sfs             = pd.DataFrame(np.round(sfs,5),columns=['f','pi','p0','pw','ps'])
		sfs['pi_nopos'] = sfs['pi'] - sfs['ps'] - sfs['pw']
		d.columns       = ['di','d0','dw','ds','di_partial','d0_partial']
	except:
		sfs             = pd.DataFrame(np.round(sfs,5),columns=['f','pi','p0','pw'])
		sfs['pi_nopos'] = sfs['pi'] - sfs['pw']
		d.columns       = ['di','d0','dw','ds']

	# Writting sfs and div file
	fname = [path + "/sfs.tsv",path + "/div.tsv"]
	sfs.to_csv(fname[0], header=True, index=False, sep="\t")
	d.to_csv(fname[1], header=True, index=False, sep="\t")
	
	idx   = [np.sort(np.random.choice(idx_all,len(daf_files),replace=True)) for i in range(1,replicas+1)]

	for i in tqdm(range(0,replicas)):
		bootstrap_sfs_and_fixations(idx[i],l_sfs,l_div,path,i)

def bootstrap_sfs_and_fixations(idx,l_sfs,l_div,path,r):

	tmp_daf          = [l_sfs[i] for i in idx]
	tmp_div          = [l_div[i] for i in idx]
	m                = np.array([50000*2000])

	sfs = np.zeros((tmp_daf[0].shape[0],tmp_daf[0].shape[1]-1))
	d = np.zeros((tmp_div[0].shape))
	for i in tmp_daf:
		sfs += i[:,1:]

	# sfs             = np.sum(tmp_daf,axis=0)[:,1:]
	f               = np.reshape(np.arange(1,sfs.shape[0]+1)/(sfs.shape[0]+1),(sfs.shape[0],1))
	sfs             = np.hstack((f,sfs))

	d               = np.sum(tmp_div,axis=0)
	d               = pd.DataFrame(d).T

	sfs             = pd.DataFrame(np.round(sfs,5),columns=['f','pi','p0','pw'])
	sfs['pi_nopos'] = sfs['pi'] - sfs['pw']
	d.columns       = ['di','d0','dw','ds']

	fname           = [path + "/sfs"+str(r)+".tsv", path + "/div"+str(r)+".tsv"]
	sfs.to_csv(fname[0], header=True, index=False, sep="\t")
	d.to_csv(fname[1], header=True, index=False, sep="\t")

def open_files(sfsFile,divFile):

	sfs = pd.read_csv(sfsFile,sep='\t').to_numpy()
	div = pd.read_csv(divFile,sep='\t').to_numpy().flatten()

	return(sfs,div)

