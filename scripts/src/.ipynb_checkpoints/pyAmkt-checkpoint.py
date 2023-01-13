from numba import njit
import numpy as np

def amkt(daf, div, xlow, xhigh):
	res = {}

	d_ratio = float(div[1] / div[0])
	# Compute alpha values and trim
	alpha = 1 - d_ratio * (daf[:,1] / daf[:,2])
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

	res['a'] = popt[0]
	res['b'] = popt[1]
	res['c'] = popt[2]

	# alpha for predicted model
	res['alpha'] = exp_model(1.0, res['a'], res['b'], res['c'])

	# Compute confidence intervals based on simulated data (MC-SOERP)
	vcov = pd.concat([pd.DataFrame([0] * 4).transpose(),
					  pd.concat([pd.DataFrame([0] * 4), pd.DataFrame(model[1])], axis=1, ignore_index=True)],
					 axis=0, ignore_index=True)
	vcov = vcov.iloc[0:4, :].values

	simpars = np.random.multivariate_normal(mean=[1.0, res['a'], res['b'], res['c']], cov=vcov, size=10000,
											check_valid='ignore')

	res['ciLow'], res['ciHigh'] = np.quantile([exp_model(x[0], x[1], x[2], x[3]) for x in simpars], [0.025, 0.975])

	return res

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
			out[i,:] = np.zeros(app.shape[0])

	return out
