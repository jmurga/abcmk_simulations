############################################
using Distributed
addprocs(20)
@everywhere include("/home/jmurga/mkt/202004/scripts/src/summaryParser.jl")
@everywhere using Analytical, CSV, DataFrames, ProgressMeter, JLD2

####
adap = Analytical.parameters(N=1000,n=500,gamNeg=-457,gL=10,gH=500,bRange=permutedims(push!(collect(0.1:0.025:0.975),0.999)),dac = [2,4,5,10,20,50,200,500,700])

convolutedSamples = Analytical.binomialDict()
Analytical.binomOp!(adap,convolutedSamples.bn);

@time df = Analytical.ratesToStats(param = adap,convolutedSamples=convolutedSamples,gH=[600,500,400,300],gL=collect(5:10),gamNeg=[-200,-300,-400,-500,-1000],iterations = 10^2,shape=adap.al,output="/home/jmurga/test.jld2");
run(`rm /home/jmurga/test.jld2`)

@time df = Analytical.ratesToStats(param = adap,convolutedSamples=convolutedSamples,gH=collect(300:600),gL=collect(5:10),gamNeg=collect(-1000:-200),iterations = 10^2,shape=adap.al,output="/home/jmurga/test.jld2");

########################################
using JLD2,Analytical,DataFrames,CSV,ProgressMeter

file = jldopen("/home/jmurga/ratesBgs.jld2","r")

adap = Analytical.parameters(N=1000,n=661,gamNeg=-457,gL=10,gH=500,bRange=permutedims(push!(collect(0.025:0.025:0.975),0.999)))
adap.dac = file[string(adap.N)*"/"*string(adap.n)*"/dac"]

analysis = ["wg" "vips" "nonvips"]

function mapDistributionData(;analysis,popSize,sample,summstatSize,replicas,h5file,output,cumulative=true)

	#=for (x,y) in zip(popSize,sample)
		print(x)=#
	adap = Analytical.parameters(N=popSize,n=sample,al=0.184)
	data = Analytical.parseSfs(sample=661,data="/home/jmurga/mkt/202004/rawData/tgp/" * analysis * ".tsv" ,dac=adap.dac)
	sCumu = Analytical.cumulativeSfs(data[4])
	divergence = data[5]

	out = "/home/jmurga/mkt/202004/rawData/summStat/tgp/"*string(analysis)
	run(`mkdir -p $out`)


	α = @. 1 - (divergence[2]/divergence[1] * sCumu[:,2]/sCumu[:,3])[adap.dac]
	α = round.(α,digits=5)
	CSV.write(out * "/alpha_"*analysis*".tsv",DataFrame(repeat(α',2)),delim='\t',header=false)

	@showprogress for i=1:replicas

		summstat, models = Analytical.summaryStatsFromRates(param = adap,rates=h5file,divergence=data[3],sfs=data[2],summstatSize=summstatSize)

		dataInfer = hcat(round.(summstat[:,1:3],digits=5),round.(Array(models[[:gamNeg,:al]]),digits=5),summstat[:,4:end])

		CSV.write(out *"/"*analysis* "_" * string(i) * ".tsv",DataFrame(dataInfer),delim='\t',header=false,append=true)

	end	
end

###################################
using Distributed
addprocs(10)
@everywhere include("/home/jmurga/mkt/202004/scripts/src/summaryParser.jl")
@everywhere using Analytical, CSV, DataFrames, ProgressMeter, JLD2, StatsBase

#=function summaryStatsFromRates(;param::parameters,rates::JLD2.JLDFile,divergence::Array,sfs::Array,summstatSize::Int64,replicas::Int64)

    tmp     = rates[string(param.N) * "/" * string(param.n)]
	idx     = StatsBase.sample.(fill(1:size(tmp["neut"],1),replicas),fill(summstatSize,replicas),replace=false)

    models  = Array.(view.(fill(tmp["models"],replicas),idx,:));
    #=alphas  = Array.(view.(fill(tmp["alphas"],replicas),idx,:));=#
    neut    = Array.(view.(fill(tmp["neut"],replicas),idx,:));
    sel     = Array.(view.(fill(tmp["sel"],replicas),idx,:));
    dsdn    = Array.(view.(fill(tmp["dsdn"],replicas),idx,:));
    #=alphas  = Array.(view.(fill(tmp["alphas"],replicas),idx,:));=#
	
	expectedValues = SharedArray{Float64,3}(summstatSize,	size(tmp["dac"],1) + 5, size(idx,1))
	#=@sync @distributed for i in eachindex(idx)
		expectedValues[:,:,i] = ratesToSummaries(alphas[i],models[i],sfs[i],divergence[i],neut[i],sel[i],dsdn[i]);
	end=#
    expectedValues =  pmap(ratesToSummaries,models,sfs,divergence,neut,sel,dsdn);

	return(expectedValues)
en
	
	re
@everwhfunction ratesToSummaries(m::Array,s::Array,d::Array,nt::Array,sl::Array,x::Array)
    ds      = x[:,1]
    dn      = x[:,2]
    dweak   = x[:,3]
    dstrong = x[:,4]
    gn      = abs.(m[:,4])
    sh      = round.(m[:,end-1],digits=5)

	alxSummStat, alphasDiv, expectedDn, expectedDs, expectedPn, expectedPs = sampledAlpha(d=d,afs=s,λdiv=[ds,dn,dweak,dstrong],λpol=[permutedims(nt),permutedims(sl)])
	expectedValues = hcat(al,gn,sh,alxSummStat)
end
=#
function mapDistribution(analysis,popSize,sample,summstatSize,replicas,h5file,simulations,output,slim=true)

	popSize=1000;sample=500;output="isolation_bgs_v1";analysis="isolation_0.4_0.1_0.2";replicas=500;summstatSize=10^5

	adap = Analytical.parameters(N=popSize,n=sample,al=0.184)
	adap.dac = h5file[string(popSize) * "/" * string(sample) * "/dac"]

	out = "/home/jmurga/mkt/202004/rawData/summStat/" * output * "/"*string(analysis)
	run(`mkdir -p $out`)

	r = collect(1:replicas)
	if slim 
		sFile = @. "/home/jmurga/mkt/202004/rawData/simulations/" * simulations *"/" * analysis * "/sfs" * string(r) * ".tsv"
		dFile = @. "/home/jmurga/mkt/202004/rawData/simulations/" * simulations *"/" * analysis * "/div" * string(r) * ".tsv"
		header = fill(false,replicas)
	else
		sFile = fill("/home/jmurga/mkt/202004/rawData/simulations"* simulations *"/" * analysis * "/sfs.tsv",replicas)
		dFile = fill("/home/jmurga/mkt/202004/rawData/simulations/" * simulations * analysis * "/div.tsv",replicas)
		header = fill(true,replicas)
	end

	tmp = openSfsDiv.(sFile,dFile,fill(adap.dac,replicas),header)

	sfs = [tmp[i][1] for i=1:replicas]
	d   = [tmp[i][2] for i=1:replicas]
	α   = [tmp[i][3] for i=1:replicas]

	@time summstat = Analytical.summaryStatsFromRates(param = adap,rates=h5file,divergence=d,sfs=sfs,summstatSize=summstatSize,replicas=replicas)

	@showprogress for i in eachindex(α)
		CSV.write(out * "/alpha" * string(i) * ".tsv",DataFrame(repeat(α[i]',2)),delim='\t',header=false)
		CSV.write(out *"/" * analysis * "_" * string(i) * ".tsv",DataFrame(summstat[i]),delim='\t',header=false)
	end
end

analysis = ["noDemog_0.4_0.1_0.2" "noDemog_0.4_0.2_0.2" "noDemog_0.4_0.3_0.2" "noDemog_0.4_0.1_0.4" "noDemog_0.4_0.2_0.4" "noDemog_0.4_0.3_0.4" "noDemog_0.4_0.1_0.8" "noDemog_0.4_0.2_0.8" "noDemog_0.4_0.3_0.8" "noDemog_0.4_0.1_0.999"  "noDemog_0.4_0.2_0.999" "noDemog_0.4_0.3_0.999"]
n = size(analysis,2)
h5file = jldopen("/home/jmurga/ratesBgs.jld2","r")
map(mapDistribution,analysis,fill(1000,n),fill(500,n),fill(10^5,n),fill(100,n),fill(h5file,n),fill("noDemog/ne500/",n),fill("noDemog_bgs_v1",n))
map(mapDistribution,analysis,fill(1000,n),fill(500,n),fill(10^5,n),fill(100,n),fill(h5file,n),fill("noDemog/ne500/",n),fill("noDemog_bgs_v1_sfs",n),fill(false,n))

######################
#=function mapDistribution(analysis,popSize,sample,summstatSize,replicas,h5file,simulations,output,cumulative=true)

	#=for (x,y) in zip(popSize,sample)
		print(x)=#
	adap = Analytical.parameters(N=popSize,n=sample,al=0.184)
	adap.dac = h5file[string(popSize) * "/" * string(sample) * "/dac"]

	sfs,divergence = readSfsDiv(analysis=analysis,path="/home/jmurga/mkt/202004/rawData/simulations/" * simulations,cumulative=cumulative)
	s = sum(sfs,dims=2)
	d = [sum(divergence[1:2])]

	out = "/home/jmurga/mkt/202004/rawData/summStat/" * output *"/"*string(analysis)
	run(`mkdir -p $out`)

	@showprogress for i=1:replicas

		summstat, models = summaryStatsFromRates_v2(param = adap,rates=h5file,divergence=d,sfs=s,summstatSize=summstatSize)

		tmp = hcat(summstat[:,1:3],Array(models[[:gamNeg,:al]]),summstat[:,4:end])
		#=tmp = tmp[tmp[:,1] .> 0 ,:]
		tmp = tmp[tmp[:,2] .> 0 ,:]
		tmp = tmp[tmp[:,3] .> 0 ,:]=#

		CSV.write(out *"/"*analysis* "_" * string(i) * ".tsv",DataFrame(tmp),delim='\t',header=false,append=false)
	end	

	α = @. 1 - (divergence[2]/divergence[1] * sfs[:,1]/sfs[:,2])[adap.dac]
	α = round.(α,digits=5)
	CSV.write(out * "/alpha_"*analysis*".tsv",DataFrame(repeat(α',2)),delim='\t',header=false)
end

function summaryStatsFromRates_v2(;param::Analytical.parameters,rates::JLD2.JLDFile,divergence::Array,sfs::Array,summstatSize::Int64)

    tmp     = rates[string(param.N) * "/" * string(param.n)]
    idx     = StatsBase.sample(1:size(tmp["neut"],1),summstatSize,replace=false)

    models  = @. round(abs(tmp["models"][idx,:]),digits=5)
    alphas  = tmp["alphas"][idx,:]
    neut    = permutedims(convert(Array,tmp["neut"]))[:,idx]
    sel     = permutedims(convert(Array,tmp["sel"]))[:,idx]
    dsdn    = permutedims(convert(Array,tmp["dsdn"]))[:,idx]
    alphas  = convert(Array,tmp["alphas"])[idx,:]
	
    ds      = dsdn[1,:]
    dn      = dsdn[2,:]
    dweak   = dsdn[3,:]
    dstrong = dsdn[4,:]

	alxSummStat, alphasDiv, expectedDn, expectedDs, expectedPn, expectedPs = Analytical.sampledAlpha(d=divergence,afs=sfs[tmp["dac"]],λdiv=[ds,dn,dweak,dstrong],λpol=[neut,sel])

	expectedValues = hcat(alphas,alxSummStat)

	return(expectedValues,models)
end=#