# from numba import njit
import numpy as np
from scipy import optimize
from scipy import stats
import pandas as pd 
from fisher import pvalue

# DAF SHOULD BE A NUMPY ARRAY CONTAINING FREQ, PI AND P0 IN THAT ORDER
# DIV SHOULD BE A NUMPY ARRAY CONTAINING DI, D0, MI, M0
# CUMULATIVESFS. DESCRIBED AT HTTPS://STATIC-CONTENT.SPRINGER.COM/ESM/ART%3A10.1038%2FS41559-019-0890-6/MEDIAOBJECTS/41559_2019_890_MOESM1_ESM.PDF USING ALL THE INFORMATION ABOVE THE FREQUENCY. THE FIRST CATEGORY AT SFS INCLUDE PI

def amkt(daf, div, xlow=0, xhigh=1):
	output = {}

	dRatio = float(div[1] / div[0])
	# Compute alpha values and trim
	alpha = 1 - dRatio * (daf[:,1] / daf[:,2])
	trim = ((daf[:,0] >= xlow) & (daf[:,0] <= xhigh))

	# Two-step model fit:
	# First bounded fit:
	try:
		popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim],bounds=([-1, -1, 1], [1, 1, 10]))
		# print('fit initial')
	except:
		# print('could not fit initial')
		popt = None
		pcov = None


	# Second fit using initially guessed values or unbounded fit:
	try:
		popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim], p0=popt,method='lm')
		# print('Fit: lm')
	except:
		try:
			popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim], p0=popt,  method='trf')
			# print('Fit: trf')
		except:
			try:
				popt, pcov = optimize.curve_fit(exp_model, daf[:,0][trim], alpha[trim], p0=popt, method='dogbox')
				# print('Fit: dogbox')
			except:
				popt=None

	if popt is None:
		output = {'a': np.nan,'b':np.nan,'c': np.nan,'alpha':np.nan,'ciLow':np.nan,'ciHigh':np.nan}

	else:
		output['a'] = popt[0]
		output['b'] = popt[1]
		output['c'] = popt[2]

		# alpha for predicted model
		output['alpha'] = exp_model(1.0, output['a'], output['b'], output['c'])
		
		# Compute confidence intervals based on simulated data (MC-SOERP)
		vcov = pd.concat([pd.DataFrame([0] * 4).transpose(),
						  pd.concat([pd.DataFrame([0] * 4), pd.DataFrame(pcov)], axis=1, ignore_index=True)],
						 axis=0, ignore_index=True)
		vcov = vcov.iloc[0:4, :].values

		try:
			simpars = np.random.multivariate_normal(mean=[1.0, output['a'], output['b'], output['c']], cov=vcov, size=10000, check_valid='raise')  # same as R implementation
			output['ciLow'], output['ciHigh'] = np.quantile([exp_model(x[0], x[1], x[2], x[3]) for x in simpars], [0.025, 0.975])
		except:
			output['ciLow'], output['ciHigh']= np.nan, np.nan

	return output,alpha[trim]

# @njit
def exp_model(f_trimmed, a, b, c):
	return a + b * np.exp(-c * f_trimmed)

# @njit
def cumulative_sfs(x):

	f       = x[:,0]
	sfsTemp = x[:,1:]
	
	out       = np.empty_like(x)
	out[0,0]  = f[0]
	out[0,1:] = np.sum(sfsTemp,axis=0)

	for i in range(1,out.shape[0]):

		app = out[i-1,1:] - sfsTemp[i-1,:]

		if np.sum(app) > 0.0:
			out[i,0]  = f[i]
			out[i,1:] = app
		else:
			out[i,0]  = f[i]
			out[i,1:] = np.zeros(app.shape[0])

	return out

def reduce_sfs(x,bins):

	bins = (bins*2) -1 
	f = x[:,0]
	sfs = x[:,1:]
	
	b = np.arange(0,1,1/bins)
	inds = np.digitize(f,b,right=True)
	out  = np.zeros((bins,x.shape[1]))
	out[:,0]  = np.unique(inds)
	
	sfsGrouped = np.hstack([np.reshape(inds,(inds.shape[0],1)),sfs])
	for i in np.unique(inds):
		out[out[:,0]==i,1:] = np.sum(sfsGrouped[sfsGrouped[:,0] == i,1:],axis=0) 
		
	out[:,0] = b
	return(out)

def imputedMK(daf, divergence,l,h=None,m=None):

	output = {}

	pi = np.sum(daf[:,1])
	p0 = np.sum(daf[:,2])
	di = divergence[0]
	d0 = divergence[1]


	toFix = 0
	deleterious = 0
	piHigh = 0
	p0High = 0
	### Estimating slightly deleterious with pi/p0 ratio
	fltLow = (daf[:, 0] <= l)
	piLow   = daf[fltLow][:,1].sum()
	p0Low   = daf[fltLow][:,2].sum()

	if (h is None):
		fltInter = (daf[:, 0] >= l) & (daf[:, 0] <= 1)
		piInter = daf[fltInter][:, 1].sum()
		p0Inter = daf[fltInter][:, 2].sum()
	else:
		fltInter = (daf[:, 0] >= l) & (daf[:, 0] < h)
		piInter = daf[fltInter][:, 1].sum()
		p0Inter = daf[fltInter][:, 2].sum()

		fltHigh = (daf[:, 0] >= h)
		piHigh  = daf[fltHigh][:, 1].sum()
		p0High  = daf[fltHigh][:, 2].sum()

		toFix = piHigh - (piInter*p0High/p0Inter)
		if toFix < 0:
			toFix = 0
		di = round(di + abs(toFix),3)


	ratioP0       = p0Low / p0Inter
	deleterious   = piLow - (piInter * ratioP0)
	piNeutral     = round(pi - deleterious - toFix,3)

	output['alpha'] = round(1 - (((pi - deleterious) / p0) * (d0 / di)),3)
	# output[0] = round(1 - ((piNeutral/p0) * (d0 / di)),3)

	if(m is not None):
		m0 = m[1];mi = m[0]
		# ## Estimation of b: weakly deleterious
		output['b'] = (deleterious / p0) * (m0 / mi)

		## Estimation of f: neutral sites
		output['f'] = (m0 * piNeutral) / (mi * p0)

		## Estimation of d, strongly deleterious sites
		output['d'] = 1 - (output['f'] + output['b'])

	oddsratio, output['pvalue'] = stats.fisher_exact([[p0, d0], [pi - deleterious, di]])

	# # divergence metrics
	# output['Ka']       = di / mi
	# output['Ks']       = d0 / m0
	# output['omega']    = output['Ka'] / output['Ks']


	## Omega A and Omega D
	# output['omegaA'] = output['omega'] * output['alpha']
	# output['omegaD'] = output['omega'] - output['omegaA']
	# r1 = round(piLow/p0Low,3)
	# r2 = round(piInter/p0Inter,3)
	# try:
	#     r3 = round(piHigh/p0High,3)
	# except:
	#     r3 = np.nan
	# return output, np.array([pi, piNeutral, p0, r1,r2,r3,di, d0])
	return output

def make_sfs_v2(data, cum=False):
	div = {'mi': data.mi.sum(),
		   'Di': data.di.sum(),
		   'm0': data.m0.sum(),
		   'D0': data.d0.sum()}

	daf = {'daf': np.arange(0.025, 1.0, 0.025),
		   'Pi': np.array(tuple(map(sum, zip(*tuple([tuple(map(int, daf.split(';'))) for daf in data.daf0f]))))),
		   'P0': np.array(tuple(map(sum, zip(*tuple([tuple(map(int, daf.split(';'))) for daf in data.daf4f])))))}

	if cum:
		daf['Pi'] = np.cumsum(daf['Pi'][::-1])[::-1]
		daf['P0'] = np.cumsum(daf['P0'][::-1])[::-1]

	return pd.DataFrame(daf, index=range(39)), pd.DataFrame(div, index=[0])

def FWW(daf, div, cutoff=0.15):
	res = {}

	P0 = daf['P0'].sum()
	Pi = daf['Pi'].sum()
	D0 = int(div['D0'])
	Di = int(div['Di'])
	m0 = int(div['m0'])
	mi = int(div['mi'])

	### Estimating alpha with Pi/P0 ratio
	PiGreater = daf[daf['daf'] > cutoff]['Pi'].sum()
	P0Greater = daf[daf['daf'] > cutoff]['P0'].sum()

	res['alpha_fww'] = 1 - ((PiGreater / P0) * (D0 / Di))
	res['pvalue_fww'] = stats.fisher_exact([[P0Greater, D0], [PiGreater, Di]])[1] 

	return res

def standardMK(sfs,divergence,m):

	output = {}

	pn = np.sum(sfs[:,1])
	ps = np.sum(sfs[:,2])
	dn = divergence[0]
	ds = divergence[1]
	
	
	output["alpha"] = round(1 - ((pn/ps) * (ds / dn)),digits=3)
	#  method = :mnnlike same results R, python two.sides
	output["pvalue"] = pvalue(ps,pn,ds,dn)

	mn = m[0]; ms = m[1]

	ka      = dn / mn
	ks       = ds / ms
	output["omega"]    = ka/ ks


	# Omega A and Omega D
	output["omegaA"] = output["omega"] * output["alpha"]
	output["omegaD"] = output["omega"] - output["omegaA"]	    

	return output
