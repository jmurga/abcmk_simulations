from numba import njit
import numpy as np
from scipy import optimize
from scipy import stats
import pandas as pd 

# DAF SHOULD BE A NUMPY ARRAY CONTAINING FREQ, PI AND P0 IN THAT ORDER
# DIV SHOULD BE A NUMPY ARRAY CONTAINING DI, D0, MI, M0
# CUMULATIVESFS. DESCRIBED AT HTTPS://STATIC-CONTENT.SPRINGER.COM/ESM/ART%3A10.1038%2FS41559-019-0890-6/MEDIAOBJECTS/41559_2019_890_MOESM1_ESM.PDF USING ALL THE INFORMATION ABOVE THE FREQUENCY. THE FIRST CATEGORY AT SFS INCLUDE PI

def amkt(daf, div, xlow, xhigh):
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
                if not popt:
                    # print('Could not fit any unbounded')
                    raise RuntimeError("Couldn't fit any method")

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

    simpars = np.random.multivariate_normal(mean=[1.0, output['a'], output['b'], output['c']], cov=vcov, size=10000, check_valid='raise')  # same as R implementation

    output['ciLow'], output['ciHigh'] = np.quantile([exp_model(x[0], x[1], x[2], x[3]) for x in simpars], [0.025, 0.975])

    return output,alpha[trim]

@njit
def exp_model(f_trimmed, a, b, c):
    return a + b * np.exp(-c * f_trimmed)

@njit
def cumulativeSfs(x):

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

@njit
def reduceSfs(x,bins):

    f = x[:,0]
    sfs = x[:,1:]
    
    b = np.arange(0,1,1/bins)
    inds = np.digitize(f,b,right=True)
    out  = np.zeros((bins,x.shape[1]))
    out[:,0]  = np.unique(inds)
    
    sfsGrouped = np.hstack([np.outputhape(inds,(inds.shape[0],1)),sfs])
    for i in np.unique(inds):
        out[out[:,0]==i,1:] = np.sum(sfsGrouped[sfsGrouped[:,0] == i,1:],axis=0) 
        
    out[:,0] = b
    return(out)

def imputedMKT(daf, div, m, cutoff=0.15):
    output = {}

    pi = np.sum(daf[:,2])
    p0 = np.sum(daf[:,1])
    di = div[1]
    d0 = div[2]
    mi = div[3]
    m0 = div[4]
    

    # divergence metrics
    output['Ka']       = di / mi
    output['Ks']       = d0 / m0
    output['omega']    = output['Ka'] / output['Ks']

    ### Estimating alpha with pi/p0 ratio
    piMinus   = daf[daf[:,0] <= cutoff][:,1].sum()
    piGreater = daf[daf[:,0] > cutoff][:,1].sum()
    p0Minus   = daf[daf[:,0] <= cutoff][:,2].sum()
    p0Greater = daf[daf[:,0] > cutoff][:,2].sum()

    ratioP0       = p0Minus / p0Greater
    deleterious   = piMinus - ( )
    piNeutral     = pi - deleterious

    output['alpha'] = 1 - (((pi - deleterious) / P0) * (d0 / di))

    ## Estimation of b: weakly deleterious
    output['b'] = (deleterious / P0) * (m0 / mi)

    ## Estimation of f: neutral sites
    output['f'] = (m0 * piNeutral) / (mi * P0)

    ## Estimation of d, strongly deleterious sites
    output['d'] = 1 - (output['f'] + output['b'])

    oddsratio, output['pvalue'] = stats.fisher_exact([[p0, d0], [pi - deleterious, di]])[1]

    ## Omega A and Omega D
    output['omegaA'] = output['omega'] * output['alpha']
    output['omegaD'] = output['omega'] - output['omegaA']

    return output
