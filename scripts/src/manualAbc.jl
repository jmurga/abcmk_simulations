using Distributed, RCall;R"""library(abc);library(data.table);library(dplyr)""";
addprocs(23)
@everywhere using Analytical, CSV, DataFrames
include("/home/jmurga/mkt/202004/scripts/src/summaryParser.jl")

adap = Analytical.parameters(N=1000,n=500,gamNeg=-457, gL=10,gH=500,Lf=2*10^5,B=0.999,alTot=0.4,alLow=0.2)
Analytical.binomOp!(adap);

runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^2,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.1_0.2",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")

runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.1_0.2",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.2_0.2",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.3_0.2",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")

runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.1_0.4",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.2_0.4",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.3_0.4",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")

runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.1_0.8",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.2_0.8",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.3_0.8",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")

runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.1_0.999",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.2_0.999",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")
runSummary(param=adap,alpha=0.5,nsamples=1,dac = [1,2,5,10,20,50,100,200,500,700,900],iterations=10^4,nthreads=5,tol=0.005,analysis="noDemog_0.4_0.3_0.999",output="/home/jmurga/mkt/202004/rawData/summStat/noDemog/",path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/")

#############################
using Distributed, RCall;
R"""library(abc);library(data.table);library(dplyr)""";
PATH="/home/jmurga/lhu/"
include("/home/jmurga/mkt/202004/scripts/src/summaryParser.jl")
addprocs(5)

@everywhere using Analytical, CSV, DataFrames
sfs = convert(Array,CSV.read(PATH * "noDemog_0.4_0.2_0.999/sfs.tsv",DataFrame))
divergence = convert(Array,CSV.read(PATH * "noDemog_0.4_0.2_0.999/div.tsv",DataFrame))

s = Analytical.cumulativeSfs(sum(sfs[:,2:3],dims=2))
d = [sum(divergence[1:2])]
adap = Analytical.parameters(N=1000,n=661,gamNeg=-457, gL=10,gH=500,Lf=2*10^5,B=0.999,alTot=0.4,alLow=0.2)
# dac = [1,2,5,10,20,50,100,200,500,700]
dac = [1,2,5,10,20,50,100,200,500]

Analytical.binomOp!(adap);
@time df = Analytical.summaryStats(param = adap,alpha = 0.5, divergence = [sum(d)],sfs = s, dac = dac,iterations = 10, scale=0.000402,shape=0.184);

# param = adap;alpha = 0.5; divergence = [sum(d)];sfs = s; dac = dac;iterations = 10;scale=0.000402;shape=0.184
sCumu = Analytical.cumulativeSfs(sfs[:,2:3])
alpha = @. 1 - (divergence[2]/divergence[1] * sCumu[:,1]/sCumu[:,2])[dac]

function ss(;param::parameters,alpha::Float64,shape::Float64=0.184,scale::Float64=0.000402,divergence::Array,sfs::Array,dac::Array{Int64,1},iterations::Int64)

	# iterations  = trunc(Int,iterations/19) + 1
	# N random prior combinations
	# fac         = rand(-2:0.05:2,iterations,2)
	# param = adap;alpha = 0.5; divergence = [sum(d)];sfs = s; dac = [1,2,4,10,20,50,100,200,500,700];iterations = 100; scale=0.000402;shape=0.184

	fac         = rand(-2:0.01:2,iterations,2)
	afac        = @. shape*(2^fac[:,1]) 
	bfac        = @. scale*(2^fac[:,2])
	alTot       = rand(collect(0.1:0.01:alpha),iterations)
	lfac        = rand(collect(0.1:0.05:0.9),iterations)
	alLow       = @. round(alTot * lfac,digits=5)
	nParam      = [param for i in 1:iterations]
	ndivergence = [divergence for i in 1:iterations]
	nSfs        = [sfs for i in 1:iterations]
	nDac        = [dac for i in 1:iterations]

	# Estimations to thread pool
	# wp  = Distributed.CachingPool(Distributed.workers())
	# tmp = Distributed.pmap(bgsIter,wp,nParam,afac,bfac,alTot,alLow,ndivergence,nSfs,nDac);
	# tmp = Distributed.pmap(bgsIter,nParam,afac,bfac,alTot,alLow,ndivergence,nSfs,nDac);
	
	out = SharedArray{Float64,3}((iterations,19,(9+size(dac,1))))
	@sync @distributed for i in eachindex(afac)
		tmp = Analytical.bgsIter(nParam[i],afac[i],bfac[i],alTot[i],alLow[i],ndivergence[i],nSfs[i],nDac[i])
		out[i,:,:] = tmp
	end

	# Output
	df = reshape(out,iterations*19,(9+size(dac,1)))
	# df  = reduce(vcat,tmp)
	# idx = vcat(1:3,3 .+ dac)
	# df  = df[:,idx]
	return df
end

@time df = ss(param = adap,alpha = 0.5, divergence = [sum(d)],sfs = s, dac = [2,4,5,10,20,50,100,200,500,700],iterations = 10^5, scale=0.000402,shape=0.184);

@time df2 = Analytical.summaryStats(param = adap,alpha = 0.5, divergence = [sum(d)],sfs = s, dac = [2,4,5,10,20,50,100,200,500,700],iterations = 10^5, scale=0.000402,shape=0.184);

