import Analytical, CSV, DataFrames
using RCall;
R"""library(ggplot2);library(data.table);library(dplyr)""";

function readSfsDiv(;analysis,path="/home/jmurga/mkt/202004/rawData/simulations/",sample=nothing)

	
	if sample != nothing
		sfs = convert(Array,CSV.read(path * analysis * "/sfs" * string(sample) * ".tsv"))
		div = convert.(Int64,convert(Array,CSV.read(path * analysis * "/div" * string(sample) * ".tsv")))
	else
		sfs = convert(Array,CSV.read(path * "/sfs.tsv",DataFrame))
		div = convert.(Int64,convert(Array,CSV.read(path * "/div.tsv",DataFrame)))
	end

	s = convert.(Int64,Analytical.cumulativeSfs(sfs[:,2:end])[:,1:2])

	return s, div
end		

function runSummary(;param,alpha,nsamples,iterations,bins,nthreads,dac,analysis,output,sample=nothing,path="/home/jmurga/mkt/202004/rawData/summStat/",parallel="/home/jmurga/.conda/envs/abcmk/bin/parallel",gzip="/usr/bin/gzip")

	# param=adap;alpha=0.5;nsamples=1;dac=[1,2,5,10,20,50,100,250,500,750,900];iterations=10^2;bins=1321;nthreads=5;analysis="noDemog_0.4_0.3_0.999";path="/home/jmurga/mkt/202004/rawData/summStat/"

	out  = output * "/" * analysis
	
	dir = `mkdir -p $out`
	run(dir)
	
	# Making alpha
	if sample !=nothing
		sfs, d = readSfsDiv(analysis=analysis,sample=sample);
	else
		sfs, d = readSfsDiv(analysis=analysis,path=path);
	end

	rSfs = Analytical.reduceSfs(sfs,bins)'
	α    = @. round(1 - (d[2]/d[1] * rSfs[:,1]/rSfs[:,2]),digits=5)[dac]'

	# Writting and gzip alpha
	CSV.write(out * "/" * "alpha_" * analysis * ".tsv", DataFrames.DataFrame(α), delim='\t', append=false,header=false)
	fAlpha = out * "/" * "alpha_" * analysis * ".tsv"
	gz = `$gzip -f $fAlpha`;
	run(gz);

	for i in 1:nsamples

		# Making summary statistics
		@time df = Analytical.summaryStats(param = adap,alpha = alpha, divergence = [sum(d)], sfs = sum(sfs,dims=2),bins = bins , dac = dac,iterations = iterations);
		
		df = df[df[:,3] .> 0,:]

		p = df[:,1:3]
		s = df[:,4:end]
		out = rcopy(R"""out = suppressWarnings(abc(param=$p,target=$α,sumstat=$s,method='rejection',tol=0.001))[['unadj.values']]""")

		# if sample !=nothing
			# CSV.write(path * "/" * "bins" * string(bins) * "/" * analysis * "/" * analysis * "_" * string(s) *  "_" * string(i) * ".tsv", DataFrames.DataFrame(out), delim='\t', append=false,header=false)
		# else
			# CSV.write(path * "/" * "bins" * string(bins) * "/" * analysis * "/" * analysis * "_" * string(i) * ".tsv", DataFrames.DataFrame(out), delim='\t', append=false,header=false)

		CSV.write(output * "/" * analysis * "/" * analysis * "_" * string(i) * ".tsv", DataFrames.DataFrame(df), delim='\t', append=false,header=false)
		CSV.write(output * "/" * analysis * "/" * analysis * "_posterior_" * string(i) * ".tsv", DataFrames.DataFrame(out), delim='\t', append=false,header=false)
	end
	
	# Gzip summaries
	fp = output * "/" * analysis * "/"
	if sample !=nothing
		f = fp .* "/" .* analysis .* "_" .* string(s) .* "_" .* string.(collect(1:nsamples)) .* ".tsv"
	else
		f = fp.* "/" .* analysis  .* "_" .* string.(collect(1:nsamples)) .* ".tsv"
		f = hcat(f,fp.* "/" .* analysis  .* "_posterior_" .* string.(collect(1:nsamples)) .* ".tsv")
	end

	gz = `$parallel -j$nthreads $gzip -f '{''}' ::: $f`;
	run(gz);
end

function discussNe(;alpha::Float64,alphaLow::Float64,bgs::Float64,pSize::Array{Int64,1},nSize::Int64,l::Int64)

	# alpha = 0.4;alphaLow = 0.2;bgs = 0.999; pSize = [100,1000,10000,100000];nSize = 100;l = 2*10^5

	out = Array{Float64}(undef,size(pSize,1),14)
	alphas  = []

	for i in eachindex(pSize)
		println(pSize[i])

		adap = Analytical.parameters(N=pSize[i],n=nSize,gam_neg=-457,Lf=l, B=bgs,gL=10,gH=500,alTot=alpha,alLow=alphaLow,bRange=[0.2 0.4 0.8 0.999])

		Analytical.binomOp!(adap)

		B = adap.B
		Analytical.set_theta_f!(adap)
		theta_f = adap.theta_f
		adap.B = 0.999
		Analytical.set_theta_f!(adap)
		Analytical.setPpos!(adap)
		adap.theta_f = theta_f
		mu=theta_f/(4*adap.N)
		adap.B = B
		x,y = Analytical.analyticalAlpha(param=adap)
		asymp = Analytical.asympFit(x)[2]


		r1 = hcat(theta_f,mu,adap.pposL,adap.pposH,adap.alLow,adap.alTot,round(x[1],digits=5),round(x[end],digits=5),round(y[1],digits=5),round(y[end],digits=5),asymp,adap.B,adap.N,adap.n)


		out[i,:] = r1

		tmp =DataFrame(alphas_nopos=y,alpha_pos=x,f=collect(1:adap.nn-1),B=adap.B)
		tmp[:Ne] = adap.N
		tmp[:Sample] = adap.n

		push!(alphas,tmp)
	end

	tmp = reduce(vcat,alphas)


	out = sort(DataFrame(out),12)
	rename!(out, [Symbol("θ"),Symbol("μ"),Symbol("pposL"),Symbol("pposH") ,Symbol("alphaW"), Symbol("alpha"),Symbol("α(x) 1"),Symbol("α(x) end"),Symbol("α(x nopos) 1"),Symbol("α(x nopos) end"),Symbol("α(asymptotic)"),Symbol("B"),Symbol("Ne"),Symbol("Sample")])
	return(out,tmp)
end
	
function readSimulations(;analysis::String,b::Float64,N::Int64,n::Int64,dac::Array{Int64,1},output::String="/home/jmurga/mkt/202004/results/simulations/alphasModel",path::String="/home/jmurga/mkt/202004/rawData/simulations/")
		# analysis="tennesen_0.4_0.3_0.999"
		# path="/home/jmurga/mkt/202004/rawData/simulations/"
		# output = "/home/jmurga/mkt/202004/results/simulations/alphasModel"
		# n = 661
		# N = 5000
		# bn = (2*n)-1
		# b=0.999
		# dac=append!([1,2,4,5,10,25,50,75],collect(100:50:900))
		# dac=[1,2,4,5,10,20,50,250,500,100]

		bn = (2*n)-1
			
		f = path * "/" * split(analysis,"_")[1] * "/" *analysis
		outPlot = output .* "/" .* analysis .* ".svg"

		tmp = parse.(Float64,split(f,"_")[2:end])
		α   = tmp[1]
		αW  = tmp[2]
		bgs = tmp[3]

		# Open files and make inputs
		sfs = convert(Array,DataFrame!(CSV.File(f * "/sfs.tsv")))
		divergence = convert(Array,DataFrame!(CSV.File(f * "/div.tsv")))

		sCumu = convert.(Int64,Analytical.cumulativeSfs(sfs[:,2:end]))
		# sCumu = sfs[:,2:end]
		d = [convert(Int64,sum(divergence[1:2]))]

		# Estimating input alpha_x
		rSfs         = Analytical.reduceSfs(sfs[:,2:end],bn)'
		rSfsCumu     = Analytical.reduceSfs(sCumu,bn)'
		alpha        = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'
		
		alphaCumu    = @. round(1 - divergence[2]/divergence[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)[dac]

		# Estimating sampled alpha_x
		br = convert(Array,append!(collect(0.1:0.05:0.95),0.999)')
		param = Analytical.parameters(N=N,n=n,B=b,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=α,alLow=αW,Lf=2*10^5,bRange = br)

		Analytical.binomOp!(param)

		j = param.B
		Analytical.set_theta_f!(param)
		theta_f = param.theta_f
		param.B = 0.999
		Analytical.set_theta_f!(param)
		Analytical.setPpos!(param)
		param.theta_f = theta_f
		param.B = j

		tmp = []
		for i in 1:100
			x2,y2,z2 = Analytical.alphaByFrequencies(param,d,sum(rSfsCumu[:,1:2],dims=2),bn,0.999,dac)
			push!(tmp,z2)
		end
		out = reduce(vcat,tmp)[:,4:end]
		t = out'

		df = R"""x=$dac;x = length(x);
			d1 = as.data.table($t);
			d1[['f']] = seq(1,x);
			d1 = melt(d1, id.vars=c('f'));
			d1 = d1 %>% group_by(f) %>% summarize(lower=quantile(value,probs=0.05),upper=quantile(value,probs=0.95),mean=mean(value));
			d1[['Estimation']] = 'Analytical';

			d2 = as.data.table($alphaCumu);
			names(d2) = 'input'
			d2[['f']] = seq(1,x)
			d2[['Estimation']] = 'Input alpha(x)'
			d2 = melt(d2, id.vars=c('f','Estimation'))
			df = merge(d1,d2,'f') %>%as.data.table

			#df = (d1,d2)
			df[['analysis']]=$analysis
			df[['B']]=$bgs
			df[['alphaW']]=$αW
			p = ggplot(df) + geom_errorbar(data=df,aes(x=f, ymin=lower, ymax=upper),color="gray", width=1) + geom_point(aes(x=f,y=value,colour='Input alpha(x)'),size=1.5) + geom_line(aes(x=f,y=value,colour='Input alpha(x)'),size=0.5) + scale_color_manual(values='#ab2710') + theme_bw();p
			df
		"""
	return(rcopy(df))
	# CSV.write("/home/jmurga/mkt/202004/results/simulations/alphasModel/alphas" * split(f,"/")[end] *".tsv", df, delim='\t';header=true)=#
end
