using Analytical, CSV, DataFrames, ProgressMeter, RCall, GZip

R"""library(ggplot2);library(data.table);library(dplyr)""";

function readSfsDiv(;analysis,path="/home/jmurga/mkt/202004/rawData/simulations/",cumulative=true,sample=nothing)

	if sample != nothing
		sample = sample * 2
		sfs = convert(Array,CSV.read(path * "/" * analysis * "/sfs" * string(sample) * ".tsv",DataFrame))
		div = convert.(Float64,convert(Array,CSV.read(path  * "/" * analysis  * "/div.tsv",DataFrame)))
	else
		sfs = convert(Array,CSV.read(path * "/" * analysis * "/sfs.tsv",DataFrame))
		div = convert.(Float64,convert(Array,CSV.read(path * "/" * analysis * "/div.tsv",DataFrame)))
	end

	if cumulative
		s = convert.(Int64,Analytical.cumulativeSfs(sfs[:,2:end])[:,1:2])
	else
		s = convert.(Int64,sfs[:,2:end])[:,1:2]
	end

	return s, div
end		

function runSummaries(;param,nsamples,iterations,nthreads,dac,tol,analysis,output,sh=[500],sl=[10],shape=0.184,sample=nothing,fixed=false,path="/home/jmurga/mkt/202004/rawData/summStat/",parallel="/home/jmurga/.conda/envs/abcmk/bin/parallel",gzip="/usr/bin/gzip")

	#=param=adap;alpha=0.2;nsamples=1;dac=[2,5,10,20,50,100,200,500,700];iterations=10^2;bins=999;nthreads=5;analysis="noDemog_0.4_0.1_0.999";output="/home/jmurga/mkt/202004/rawData/summStat/noDemogB";path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/";sample=nothing=#
	
	out  = output * "/" * analysis
	pattern = split(analysis,"_")[1]

	dir = `mkdir -p $out`
	run(dir)
	
	# Making alpha
	if sample !=nothing
		sfs, d = readSfsDiv(analysis=analysis,path=path,sample=sample);
	else
		sfs, d = readSfsDiv(analysis=analysis,path=path);
	end

	α    = @. round(1 - (d[2]/d[1] * sfs[:,1]/sfs[:,2]),digits=5)[dac]'
	α = repeat(α,2)
	al = α[1,:] 
	amk = round(Analytical.asympFit(α[1,:])[1],digits=2)
	println(amk)
	
	# Writting and gzip alpha
	CSV.write(out * "/" * "alpha_" * analysis * ".tsv", DataFrames.DataFrame(α), delim='\t', append=false,header=false)
	fAlpha = out * "/" * "alpha_" * analysis * ".tsv"
	
	for i in nsamples[1]:nsamples[2]

		# Making summary statistics
		@time df = Analytical.summaryStats(param = param,amk = amk, divergence = [sum(d)], sfs = sum(sfs,dims=2),gH=sh,gL=sl,dac = dac,iterations = iterations,shape=param.al,scale=abs(param.al/param.gamNeg),fixed=fixed);

		df = rcopy(R"""na.omit($df)""")
		p = df[:,1:3]
		s = df[:,4:end]

		outNN = rcopy(R"""out = suppressWarnings(abc(param=$p,target=$al,sumstat=$s,method='ridge',tol=$tol))[1:2]""")

		CSV.write(output * "/" * analysis * "/" * analysis * "_" * string(i) * ".tsv", DataFrames.DataFrame(df), delim='\t', append=false,header=false)
		CSV.write(output * "/" * analysis * "/rejection_" * string(i) * ".tsv", DataFrames.DataFrame(outNN[:unadj_values]), delim='\t', append=false,header=false)
		CSV.write(output * "/" * analysis * "/neuralnet_" * string(i) * ".tsv", DataFrames.DataFrame(outNN[:adj_values]), delim='\t', append=false,header=false)

		ss = length(dac)
		posterior = out *"/" * analysis * "_" * string(i) * ".tsv";
		b = out * "/out_" * string(i)
		reg = `/home/jmurga/ABCreg/src/reg -p $posterior -d $fAlpha -P 3 -S $ss -t $tol -b $b`
		regLog = `/home/jmurga/ABCreg/src/reg -p $posterior -d $fAlpha -P 3 -S $ss -L -t $tol -b $b`
		run(reg);
		run(regLog);
	end
	
	# Gzip summaries
	fp = output * "/" * analysis 
	searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

	tmp = searchdir(fp,"tsv")
	f = tmp[.!occursin.("gz",tmp)]
	f = out .* "/" .* f

	gz = `$parallel -j$nthreads $gzip --force '{''}' ::: $f`;
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
	
function readSimulations(;analysis::String,N::Int64,n::Int64,dac::Array,br::Array=permutedims(push!(collect(0.1:0.05:0.95),0.999)),output::String="/home/jmurga/mkt/202004/results/simulations/model",path::String="/home/jmurga/mkt/202004/rawData/simulations/")
	# analysis="isolation_0.4_0.1_0.999"
	# path="/home/jmurga/mkt/202004/rawData/simulations/isolation/"
	# output = "/home/jmurga/mkt/202004/results/simulations/alphasModel"
	# n = 500
	# N = 5000
	# bn = (2*n)-1
	# b=0.999
	# dac=append!([1,2,4,5,10,25,50,75],collect(100:50:900))
	# dac=[5,10,20,50,100,200,500,900]

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

	param = Analytical.parameters(N=N,n=n,gamNeg=-457,gL=10,gH=500,al=0.184,be=0.000402,alTot=α,alLow=αW,Lf=2*10^5,bRange = br,dac=dac)
	convoluted = Analytical.binomialDict()	
	Analytical.binomOp!(param,convoluted.bn)

	s = convert.(Int64,Analytical.cumulativeSfs(sfs[:,2:end]))
	sfsInput = sum(s[:,1:2],dims=2)
	sCumu = s[param.dac,:]
	d = [convert(Int64,sum(divergence[1:2]))]

	# Estimating input alpha_x
	alpha        = @. round(1 - divergence[2]/divergence[1] * sCumu[:,1]/sCumu[:,2],digits=5)
	
	# Estimating sampled alpha_x
	if (size(br,1) == 1 )
		br = permutedims(br)
	end


	tablePlot = []
	@showprogress  for i in param.bRange
		# println(i)
		param.B = i
		j = param.B
		Analytical.setThetaF!(param)
		theta_f = param.thetaF
		param.B = 0.999
		Analytical.setThetaF!(param)
		Analytical.setPpos!(param)
		param.thetaF = theta_f
		param.B = j

		tmp = []
		for n in 1:1000
			x2,y2,z2 = Analytical.alphaByFrequencies(param,convoluted,d,sfsInput)
			push!(tmp,z2)
		end
		out = reduce(vcat,tmp)[:,4:end]
		t = out'
		
		nn = param.dac
		df = R"""d1 = as.data.table($t);
			d1[['f']] = $nn
			d1 = melt(d1, id.vars=c('f'));
			d1 = suppressMessages(d1 %>% group_by(f) %>% summarize(lower=quantile(value,probs=0.05),upper=quantile(value,probs=0.95),mean=mean(value)));
			d1[['Estimation']] = 'Analytical';
			d2 = as.data.table($alpha);
			names(d2) = 'input'
			d2[['f']] = $nn
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

	# p = rcopy(R"""p = ggplot($tablePlot) + geom_errorbar(data=$tablePlot,aes(x=f, ymin=lower, ymax=upper),color="gray", width=1) + geom_point(aes(x=f,y=value,colour='Input alpha(x)'),size=1.5) + geom_line(aes(x=f,y=value,colour='Input alpha(x)'),size=0.5) + scale_color_manual(values='#ab2710') + theme_bw() + facet_wrap(~B,ncol=6);p""")
	#ggsave(p,filename=$outPlot)"""

	return tablePlot
	# CSV.write("/home/jmurga/mkt/202004/results/simulations/alphasModel/alphas" * split(f,"/")[end] *".tsv", df, delim='\t';header=true)
end

function plotModel(;df,dac,title,gridColumns,output)

	p = R"""df=$df %>% as.data.table
	
	df = df[f %in% $dac]
	
	if(length($dac) < 50){
		df$f = as.factor(df$f)
	}

	df$f2 = rep(seq(1,length($dac)),length(unique(df$B)))

	p = ggplot(df) + geom_errorbar(data=df,aes(x=f, ymin=lower, ymax=upper),color='gray', width=2) + geom_point(aes(x=f,y=value,colour='Input alpha(x)'),size=1.5) + geom_line(aes(x=f2,y=value,colour='Input alpha(x)'),size=0.5) + scale_color_manual(values='#ab2710') + theme_bw() + facet_wrap(~B,ncol=$gridColumns) +  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5)) + ggtitle($title)
	
	ggsave(p,filename=$output)
	p
	"""
	return(p)
end

function abcSummaries(analysis::String="<folder>",ne::Int64=1000, samples::Int64=500,rates::JLD2.JLDFile="rates.jld2",summstatSize::Int64=100000,replicas::Int64=100,output::String="<output>",bootstrap::Bool=false,simulations::String="/home/jmurga/mkt/202004/rawData/simulations")
		
    adap     = Analytical.parameters(N=ne,n=samples)
    adap.dac = h5file[string(ne) * "/" * string(samples) * "/dac"]

    out = output * "/" * analysis
	run(`mkdir -p $out`)

	r = collect(1:replicas)
	if bootstrap
		sFile  = fill(simulations * "/" * analysis * "/sfs.tsv",replicas)
		dFile  = fill(simulations * "/" * analysis * "/div.tsv",replicas)
		b      = fill(true,replicas)
	else
		sFile  = @. simulations * "/" * split(analysis,"_")[1] * "/" * analysis * "/sfs" * string(r) * ".tsv"
		dFile  = @. simulations * "/" * split(analysis,"_")[1] * "/" * analysis * "/div" * string(r) * ".tsv"
		b      = fill(false,replicas)
	end
	
	tmp = openSfsDiv.(sFile,dFile,fill(adap.dac,replicas),b)

	sfs = [tmp[i][1] for i=1:replicas]
	d   = [tmp[i][2] for i=1:replicas]
	α   = [tmp[i][3] for i=1:replicas]

	summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,divergence=d,sfs=sfs,summstatSize=summstatSize,replicas=replicas)

	@showprogress for i in eachindex(α)
		CSV.write(out * "/alpha_" * string(i) * ".tsv",DataFrame(repeat(α[i]',2)),delim='\t',header=false)
		CSV.write(out * "/" * "summstat_" * string(i) * ".tsv",DataFrame(summstat[i]),delim='\t',header=false)
	end
end

function abcInference(analysis::String="analysis",replicas::Int64=100,P::Int64=5,S::Int64=9,tol::Float64=0.01,workers::Int64=1,parallel::Bool=true)
	
	reg = chomp(read(`which reg`,String))

	aFile = @. analysis * "/alpha_" * string(1:replicas) * ".tsv"
	sumFile = @. analysis * "/summstat_" * string(1:replicas) * ".tsv"
	out = @. analysis * "/out_" * string(1:replicas)
	if parallel
		run(`parallel -j$workers --link $reg -d "{1}" -p "{2}" -P $P -S $S -t $tol -L -b "{3}" ::: $aFile ::: $sumFile ::: $out`)
	else
		r(a,s,o) = run(`reg -d $a -p $s -P $P -S $S -t $tol -L -b $o`)
		r.(aFile,bFile,out)
	end	
end

function openSfsDiv(x,y,dac,bootstrap=false)

	sfs = CSV.read(x,DataFrame) |> Array

	if bootstrap
		sfs[:,2:end] = PoissonRandom.pois_rand.(sfs[:,2:end])
	end

	scumu = Analytical.cumulativeSfs(sfs)
	s = sum(scumu[:,2:3],dims=2)[dac]
	divergence = CSV.read(y,DataFrame) |> Array
	d = [sum(divergence[1:2])]
	α = @. round(1 - (divergence[2]/divergence[1] * scumu[:,2]/scumu[:,3])[dac],digits=5)
	return(s,d,α)	
end

function openAbcDistributions(file)
	df         = CSV.read(GZip.open(file),DataFrame,header=false)
	neg        = df[:,4:5]
	al         = df[:,1:3]

	#=mn         = mapcols(c -> mean(c),al)=#
	mp         = rcopy(R"""suppressWarnings(apply($al,2,getmap))""")
	mneg       = rcopy(R"""suppressWarnings(apply($neg,2,getmap))""")

	tmp = vcat(mp,mneg)'
	out = DataFrame(tmp)
	rename!(out,[:aw,:as,:a,:gamNeg,:shape])
	return(out)
end

function plotMap(case,replicas,path="/home/jmurga/mkt/202004/rawData/simulations/")


	for n in case
		println(n)
		analysis = split(n,"/")[end]

		tmp = parse.(Float64,split(analysis,"_")[end-2:end])
		a = tmp[1]
		w = tmp[2]
		b = tmp[3]

		out = filter(x -> occursin("post",x), readdir(n,join=true))
		out = filter(x -> !occursin(".1.",x),out)
		iterMap = vcat(openAbcDistributions.(out)...);
		insertcols!(iterMap,:analysis => fill(analysis,size(iterMap,1)),:b => fill(b,size(iterMap,1)),:w => fill(w,size(iterMap,1)),:method => fill("abc",size(iterMap,1)))

		dFiles = @. path * "/" * split(analysis,"_")[1] * "/" * analysis * "/div" * string(1:replicas) * ".tsv"
		function openTrues(d)
			tmp = Array(CSV.read(d,DataFrame))
			return(round.([tmp[3]/tmp[1],tmp[4]/tmp[1],(tmp[3]+tmp[4])/tmp[1]],digits=5))
		end

		trueAlphas = DataFrame(hcat(openTrues.(dFiles)...)',[:aw,:as,:a])


		insertcols!(trueAlphas,:analysis => fill(analysis,size(trueAlphas,1)),:b => fill(b,size(trueAlphas,1)),:w => fill(w,size(trueAlphas,1)),:method => fill("slim",size(trueAlphas,1)))

		push!(dmap,iterMap)
		push!(slimAlphas,trueAlphas)
	end

	dmap = vcat(dmap...)
	slimAlphas = vcat(slimAlphas...)
	
	df  = vcat(dmap[:,[1,2,3,6,7,8,9]],slimAlphas)
	p = rcopy(R"""df = $df %>% arrange(b,w)
		df[["analysis"]] = factor(df[["analysis"]],levels=unique(df[["analysis"]]))
		d = melt(df,id.vars=c("analysis","b","w","method"))
		p   = ggplot(d) + geom_boxplot(aes(x=variable,y=value,fill=method),alpha=0.5) + facet_wrap(~analysis,ncol=3) + fillSanMiguel() + colourSanMiguel();ggsave()""")
	return(df)
end