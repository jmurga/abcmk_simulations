using Analytical, CSV, DataFrames 

function analyticalApproach(param,cnv)

	B = param.B
	#=Analytical.setThetaF!(param)
	theta_f = param.thetaF
	param.B = 0.999
	Analytical.setThetaF!(param)
	Analytical.setPpos!(param)
	param.thetaF = theta_f
	param.B = B=#
	x,y = Analytical.analyticalAlpha(param=param,convolutedSamples=cnv)

	mu = param.thetaF/(4*param.N)

	return hcat(mu,param.pposL,param.pposH,param.alLow,param.alTot,round(y[end],digits=5),param.B)
end

function simTable(alphas,bgsValues,pSize,nSize,l)

	out = zeros(3,7,size(bgsValues,2))
	it = 1

	adap = Analytical.parameters(N=pSize,n=nSize,gamNeg=-457,Lf=l, B=0.999,gL=10,gH=500,alTot=0.4,alLow=0.2,bRange=bgsValues)
	cnv = Analytical.binomialDict()
	Analytical.binomOp!(adap,cnv.bn)

	for a in 1:size(alphas,1)
		for b in 1:size(bgsValues,2)

			adap.B = bgsValues[b]

			adap.alLow = alphas[a]*0.25
			r1 = analyticalApproach(adap,cnv)
			adap.alLow = alphas[a]*0.5
			r2 = analyticalApproach(adap,cnv)
			adap.alLow = alphas[a]*0.75
			r3 = analyticalApproach(adap,cnv)
		
			out[:,:,b] = vcat(r1,r2,r3)
	
	  end
	end
	out = vcat(eachslice(out,dims=3)...);
	out = sort(DataFrame(out),5)
	out[:,4] = round.(out[:,4],digits=2)

	rename!(out, [Symbol("bgsThetaF"),Symbol("pposL"),Symbol("pposH") ,Symbol("alphaW"), Symbol("alpha"),Symbol("estimation"),Symbol("B")])
	
	return out
end

alphas = [0.4]
bgs   = [0.2 0.4 0.8 0.999]
bgs   = [0.999] |> permutedims

simulations = simTable(alphas,bgs,parse(Int,ARGS[1]),parse(Int,ARGS[2]),parse(Int,ARGS[3]))

CSV.write("/home/jmurga/mkt/202004/rawData/simulations/" * ARGS[4] * ".tsv", simulations, delim='\t')