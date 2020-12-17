import Analytical, CSV, DataFrames
using RCall;
R"""library(ggplot2);library(data.table);library(dplyr)""";

function readSfsDiv(;analysis,path="/home/jmurga/mkt/202004/rawData/simulations/",sample=nothing)

	if sample != nothing
		sample = sample * 2
		sfs = convert(Array,CSV.read(path * "/" * analysis * "/sfs" * string(sample) * ".tsv"))
		div = convert.(Int64,convert(Array,CSV.read(path  * "/" * analysis  * "/div.tsv")))
	else
		sfs = convert(Array,CSV.read(path * "/" * analysis * "/sfs.tsv",DataFrame))
		div = convert.(Int64,convert(Array,CSV.read(path * "/" * analysis * "/div.tsv",DataFrame)))
	end

	s = convert.(Int64,Analytical.cumulativeSfs(sfs[:,2:end])[:,1:2])

	return s, div
end		

function runSummary(;param,alpha,nsamples,iterations,nthreads,dac,tol,analysis,output,sample=nothing,path="/home/jmurga/mkt/202004/rawData/summStat/",parallel="/home/jmurga/.conda/envs/abcmk/bin/parallel",gzip="/usr/bin/gzip")

	# param=adap;alpha=0.5;nsamples=1;dac=[1,2,4,5,10,20,50,150,350,700];iterations=10^2;bins=999;nthreads=5;analysis="tennesen_0.4_0.2_0.999";output="/home/jmurga/mkt/202004/rawData/summStat/tennesen";path="/home/jmurga/mkt/202004/rawData/simulations/tennesen/rescaled/"

	out  = output * "/" * analysis
	
	dir = `mkdir -p $out`
	run(dir)
	
	# Making alpha
	if sample !=nothing
		sfs, d = readSfsDiv(analysis=analysis,path=path,sample=sample);
	else
		sfs, d = readSfsDiv(analysis=analysis,path=path);
	end

	α    = @. round(1 - (d[2]/d[1] * sfs[:,1]/sfs[:,2]),digits=5)[dac]'

	# Writting and gzip alpha
	CSV.write(out * "/" * "alpha_" * analysis * ".tsv", DataFrames.DataFrame(α), delim='\t', append=false,header=false)
	fAlpha = out * "/" * "alpha_" * analysis * ".tsv"
	gz = `$gzip -f $fAlpha`;
	run(gz);

	for i in 1:nsamples

		# Making summary statistics
		@time df = Analytical.summaryStats(param = adap,alpha = alpha, divergence = [sum(d)], sfs = sum(sfs,dims=2),dac = dac,iterations = iterations);
		
		df = df[df[:,3] .> 0,:]

		p = df[:,1:3]
		s = df[:,4:end]
		outABC = rcopy(R"""out = suppressWarnings(abc(param=$p,target=$α,sumstat=$s,method='rejection',tol=$tol))[['unadj.values']]""")

		if sample !=nothing
			CSV.write(output * "/" * analysis * "/" * analysis * "_" * string(i) * "_sample_" * string(sample) * ".tsv", DataFrames.DataFrame(df), delim='\t', append=false,header=false)
			CSV.write(output * "/" * analysis * "/" * analysis * "_posterior_" * string(i) * "_sample_" * string(sample) * ".tsv", DataFrames.DataFrame(outABC), delim='\t', append=false,header=false)
		else
			CSV.write(output * "/" * analysis * "/" * analysis * "_" * string(i) * ".tsv", DataFrames.DataFrame(df), delim='\t', append=false,header=false)
			CSV.write(output * "/" * analysis * "/" * analysis * "_posterior_" * string(i) * ".tsv", DataFrames.DataFrame(outABC), delim='\t', append=false,header=false)
		end
	end
	
	# Gzip summaries
	fp = output * analysis 
	searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
	if sample !=nothing
		tmp = out .* "/" .*  searchdir(fp,"sample")
		f = tmp[.!occursin.("gz",tmp)]
	else
		tmp = searchdir(fp,"noDemog")
		tmp = tmp[.!occursin.("sample",tmp)]
		tmp = tmp[.!occursin.("alpha",tmp)]
		f = tmp[.!occursin.("gz",tmp)]
		f = out .* "/" .* f
	end

	gz = `$parallel -j$nthreads $gzip --force '{''}' ::: $f`;
	# gz = `$gzip $f`;
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
	
function readSimulations(;analysis::String,N::Int64,n::Int64,output::String="/home/jmurga/mkt/202004/results/simulations/model",path::String="/home/jmurga/mkt/202004/rawData/simulations/")
	# analysis="tennesen_0.4_0.1_0.2"
	# path="/home/jmurga/mkt/202004/rawData/simulations/"
	# output = "/home/jmurga/mkt/202004/results/simulations/alphasModel"
	# n = 500
	# N = 5000
	# bn = (2*n)-1
	# b=0.999
	# dac=append!([1,2,4,5,10,25,50,75],collect(100:50:900))
	# dac=[1,2,4,5,10,20,50,250,500,700]

	bn = (2*n)-1
		
	f = path * "/" * analysis
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
	rSfs         = sfs
	rSfsCumu     = sCumu
	alpha        = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'
	
	alphaCumu    = @. round(1 - divergence[2]/divergence[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)

	# Estimating sampled alpha_x
	br = convert(Array,append!(collect(0.1:0.05:0.95),0.999)')
	param = Analytical.parameters(N=N,n=n,gam_neg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=α,alLow=αW,Lf=2*10^5,bRange = br)
	Analytical.binomOp!(param)

	tablePlot = []
	for i in param.bRange
		# println(i)
		param.B = i
		j = param.B
		Analytical.set_theta_f!(param)
		theta_f = param.theta_f
		param.B = 0.999
		Analytical.set_theta_f!(param)
		Analytical.setPpos!(param)
		param.theta_f = theta_f
		param.B = j

		tmp = []
		for n in 1:100
			x2,y2,z2 = Analytical.alphaByFrequencies(param,d,sum(rSfsCumu[:,1:2],dims=2),bn,0.999,collect(1:param.nn-1))
			push!(tmp,z2)
		end
		out = reduce(vcat,tmp)[:,4:end]
		t = out'
		
		nn = param.nn-1
		df = R"""x=seq(1,$nn);x = length(x);
			d1 = as.data.table($t);
			d1[['f']] = seq(1,x);
			d1 = melt(d1, id.vars=c('f'));
			d1 = suppressMessages(d1 %>% group_by(f) %>% summarize(lower=quantile(value,probs=0.05),upper=quantile(value,probs=0.95),mean=mean(value)));
			d1[['Estimation']] = 'Analytical';
			d2 = as.data.table($alphaCumu);
			names(d2) = 'input'
			d2[['f']] = seq(1,x)
			d2[['Estimation']] = 'Input alpha(x)'
			d2 = melt(d2, id.vars=c('f','Estimation'))
			df = merge(d1,d2,'f') %>%as.data.table
			#df = (d1,d2)
			df[['analysis']]=$analysis
			df[['B']]=$i
			df[['alphaW']]=$αW
			# p = ggplot(df) + geom_errorbar(data=df,aes(x=f, ymin=lower, ymax=upper),color="gray", width=1) + geom_point(aes(x=f,y=value,colour='Input alpha(x)'),size=1.5) + geom_line(aes(x=f,y=value,colour='Input alpha(x)'),size=0.5) + scale_color_manual(values='#ab2710') + theme_bw();p
			df"""
		push!(tablePlot,rcopy(df))
	end
	tablePlot = reduce(vcat,tablePlot)

	# R"""p = ggplot($tablePlot) + geom_errorbar(data=$tablePlot,aes(x=f, ymin=lower, ymax=upper),color="gray", width=1) + geom_point(aes(x=f,y=value,colour='Input alpha(x)'),size=1.5) + geom_line(aes(x=f,y=value,colour='Input alpha(x)'),size=0.5) + scale_color_manual(values='#ab2710') + theme_bw() + facet_wrap(~B,ncol=6);ggsave(p,filename=$outPlot)"""

	return tablePlot
	# CSV.write("/home/jmurga/mkt/202004/results/simulations/alphasModel/alphas" * split(f,"/")[end] *".tsv", df, delim='\t';header=true)
end

function plotModel(;df,dac,title,gridColumns,output)

	p = R"""df=$df %>% as.data.table
	
	df = df[f %in% $dac]
	
	df$f = as.factor(df$f)
	
	df$f2 = rep(seq(1,length($dac)),length(unique(df$B)))

	p = ggplot(df) + geom_errorbar(data=df,aes(x=f, ymin=lower, ymax=upper),color='gray', width=2) + geom_point(aes(x=f,y=value,colour='Input alpha(x)'),size=1.5) + geom_line(aes(x=f2,y=value,colour='Input alpha(x)'),size=0.5) + scale_color_manual(values='#ab2710') + theme_bw() + facet_wrap(~B,ncol=$gridColumns) +  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5)) + ggtitle($title)
	
	ggsave(p,filename=$output)
	p
	"""
	return(p)
end