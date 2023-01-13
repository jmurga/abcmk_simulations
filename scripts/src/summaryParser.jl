using Analytical, CSV, DataFrames, ProgressMeter, RCall, GZip, StatsBase, JLD2

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
	
function modeledAlpha(;analysis::String,n::Int64,dac::Array,output::String,path::String,title::String,br::Array=permutedims(push!(collect(0.05:0.05:0.95),0.999)),randomize::Bool=false)

	N  = 1000
	
	f = path * "/" * analysis

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

	s = Analytical.cumulativeSfs(sfs)
	sfsInput = sum(s[:,2:3],dims=2)[param.dac,:]
	d = [convert(Int64,sum(divergence[1:2]))]

	# Estimating input alpha_x
	empircalAlpha = @. round(1 - divergence[2]/divergence[1] * s[:,2]/s[:,3],digits=5)[param.dac]
	
	iterations=1000
	#=out    = Array{Float64,3}(size(param.bRange,2),(size(param.dac,1) *2) + 12,iterations)=#
	if randomize
		alpha     = rand(0:0.1:0.9)
		alphaWeak = rand(0:0.05:0.9) * alpha
		bval      = rand(push!(collect(0.05:0.05:0.9),0.999))
		al        = round(0.184 * 2^rand(-2:0.05:2),digits=3)
		gH        = rand(200:1000)
		gL        = rand(5:10)
		gamNeg    = rand(-2000:-200)
		outPlot1  = output .* "/modeled_" .* analysis .* "_Random.svg"
		outPlot2  = output .* "/modeled_" .* analysis .* "_Random.jpg"
	else
		alpha     = α
		alphaWeak = αW
		bval      = bgs
		al        = 0.184
		gH        = 500
		gL        = 10
		gamNeg    = -457
		outPlot1 = output .* "/modeled_" .* analysis .* ".svg"
		outPlot2 = output .* "/modeled_" .* analysis .* ".jpg"
	end

	out    = zeros(size(param.bRange,2),(size(param.dac,1) *2) + 12,iterations)
	for i in 1:iterations
		# Each iteration solve 1 model accounting all B value in param.bRange
		@inbounds out[:,:,i] = Analytical.iterRates(fill(param,iterations)[i], fill(convoluted,iterations)[i], fill(alpha,iterations)[i], fill(alphaWeak,iterations)[i], fill(gH,iterations)[i], fill(gL,iterations)[i], fill(gamNeg,iterations)[i],fill(al,iterations)[i]);
	end

	# Reducing array
	df = vcat(eachslice(out,dims=3)...);
	
	# Saving models and rates
	models = [df[:,1:8]]
	neut   = [df[:,9:(8+size(param.dac,1))]]
	sel    = [df[:,(9+size(param.dac,1)):(8+size(param.dac,1)*2)]]
	dsdn   = [df[:,(end-3):end]]

	expectedValues =  progress_map(Analytical.samplingFromRates,models,[sfsInput],[d],neut,sel,dsdn;progress=Progress(iterations,desc="Estimating summaries "));
	tablePlot = DataFrame(hcat(df[:,1],expectedValues[1][:,6:end]))

	df = R"""d1 = as.data.table($tablePlot);

		d1 = melt(d1, id.vars=c('x1'));
		d1 = d1 %>% group_by(x1,variable) %>% summarize(upper=quantile(value,probs=0.05),lower=quantile(value,probs=0.95),mean=mean(value))
		levels(d1[['variable']]) = $dac
		d1[['Estimation']]= 'Analytical'
		names(d1)[1:2] = c('B','f')
		d2 = as.data.table($empircalAlpha);
		names(d2) = 'input'
		d2[['f']] = factor($dac)
		d2[['Estimation']] = 'Input alpha(x)'
		d2 = melt(d2, id.vars=c('f','Estimation'))
		df = merge(d1,d2,'f') %>%as.data.table
		df[['analysis']] = $analysis
		p = ggplot(df) + geom_errorbar(data=df,aes(x=f, ymin=lower, ymax=upper),color='gray', width=2) + geom_point(aes(x=f,y=value,colour='Input alpha(x)'),size=1.5) + geom_line(aes(x=as.numeric(f),y=value,colour='Input alpha(x)'),size=0.5) + scale_color_manual(values='#ab2710') + theme_bw() + facet_wrap(~B,ncol=5) +  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=0.5)) + ggtitle($title)
		ggsave(p,filename=$outPlot1,dpi=600)
		ggsave(p,filename=$outPlot2,dpi=600)
		print(p)"""
	@rget df
	return df
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
	adap.dac = rates[string(ne) * "/" * string(samples) * "/dac"]

	out = output * "/" * analysis
	run(`mkdir -p $out`)

	r = collect(1:replicas)
	if bootstrap
		sFile  = fill(simulations * "/" * split(analysis,"_")[1] *"/"* analysis * "/sfs.tsv",replicas)
		dFile  = fill(simulations * "/" * split(analysis,"_")[1] *"/"* analysis * "/div.tsv",replicas)
		b      = fill(true,replicas)
	else
		if occursin("twoepochs", analysis)
			fName =  @. simulations * "/"* analysis
		else
			fName = @. simulations * "/" * split(analysis,"_")[1] *"/"* analysis
		end
		sFile  = @. fName * "/sfs"*string(r)*".tsv"
		dFile  = @. fName * "/div"*string(r)*".tsv"
		b      = fill(true,replicas)
	end
	
	tmpOpening = openSfsDivManual.(sFile,dFile,fill(adap.dac,replicas));
	sfs = [tmpOpening[i][1] for i in eachindex(tmpOpening)]
	d = [tmpOpening[i][2] for i in eachindex(tmpOpening)]
	α = [tmpOpening[i][3] for i in eachindex(tmpOpening)]

	#Open rates
	tmp    = rates[string(adap.N) * "/" * string(adap.n)]

	#Subset index
	idx    = StatsBase.sample.(fill(1:size(tmp["neut"],1),replicas),fill(summstatSize,replicas),replace=false)

	models = Array.(map(x -> view(tmp["models"],x,:), idx));
	neut   = Array.(map(x -> view(tmp["neut"],x,:), idx));
	sel    = Array.(map(x -> view(tmp["sel"],x,:), idx));
	neut   = Array.(map(x -> view(tmp["neut"],x,:), idx));
	dsdn   = Array.(map(x -> view(tmp["dsdn"],x,:), idx));
	
	#Making summaries
	expectedValues =  progress_pmap(Analytical.samplingFromRates,models,sfs,d,neut,sel,dsdn;progress=Progress(replicas,desc="Estimating summaries "));

	w(x,name) = CSV.write(name,DataFrame(x),delim='\t',header=false);
	
	#Writting ABCreg input
	progress_pmap(w,repeat.(α,2), @. out * "/alpha_" * string(1:replicas) * ".tsv";progress= Progress(replicas,desc="Writting alphas "));
	progress_pmap(w, expectedValues, @. out * "/summstat_" * string(1:replicas) * ".tsv";progress= Progress(replicas,desc="Writting summaries "));

end

function abcInference(analysis::String="analysis",replicas::Int64=100,P::Int64=5,S::Int64=9,tol::Float64=0.01,workers::Int64=1,parallel::Bool=true)
	
	reg = chomp(read(`which reg`,String))

	aFile = @. analysis * "/alpha_" * string(1:replicas) * ".tsv"
	sumFile = @. analysis * "/summstat_" * string(1:replicas) * ".tsv"
	out = @. analysis * "/out_" * string(1:replicas)
	if parallel
		run(`parallel -j$workers --link $reg -d "{1}" -p "{2}" -P $P -S $S -t $tol -b "{3}" ::: $aFile ::: $sumFile ::: $out`)
	else
		r(a,s,o) = run(`reg -d $a -p $s -P $P -S $S -t $tol -L -b $o`)
		r.(aFile,bFile,out)
	end	
end

function openSfsDivManual(x,y,dac,bootstrap=false)

	sfs = CSV.read(x,DataFrame) |> Array

	if bootstrap
		sfs[:,2:end] = PoissonRandom.pois_rand.(sfs[:,2:end])
	end

	#=scumu = Analytical.cumulativeSfs(Analytical.reduceSfs(sfs,40))=#
	scumu = Analytical.cumulativeSfs(sfs)
	s = sum(scumu[:,2:3],dims=2)[dac]
	divergence = CSV.read(y,DataFrame) |> Array
	d = [sum(divergence[1:2])]
	α = @. round(1 - (divergence[2]/divergence[1] * scumu[:,2]/scumu[:,3])[dac],digits=5)
	return(s,d,permutedims(α))
end

function openAbcDistributions(file)
	df         = CSV.read(GZip.open(file),DataFrame,header=false)
	neg        = df[:,4:5]
	al         = df[:,1:3]

	#=mn         = mapcols(c -> mean(c),al)=#

	mp         = rcopy(R"""suppressWarnings(apply($al,2,getmap))""")

	cnames = [:aw,:as,:a,:gamNeg,:shape]
	mneg       = rcopy(R"""suppressWarnings(apply($neg,2,getmap))""")

	tmp = vcat(mp,mneg)'
	out = DataFrame(tmp)
	rename!(out,cnames)
	return(out)
end

function plotMapSimulations(case,replicas,output,title,simulations="/home/jmurga/mkt/202004/rawData/simulations/")

	slimAlphas = []
	dmap       = []

	for n in case
		println(n)
		analysis = split(n,"/")[end]

		if(occursin("twoepochs",analysis))
			a  = 0.4
			aw = 0.1
			b  = 0.999
			if occursin("eq",analysis)
				fName = "/home/jmurga/mkt/202004/rawData/simulations/twoepochs_exp/twoepochs_eq"
			else
				fName = @. simulations * "/" * analysis
			end
			facetCols = 2
			@rput facetCols
			
		else
			tmp = parse.(Float64,split(analysis,"_")[end-2:end])
			a   = tmp[1]
			aw  = tmp[2]
			b   = tmp[3]
			fName = @. simulations * "/" * split(analysis,"_")[1] * "/" * analysis
			facetCols =3
			@rput facetCols
		end

		out = filter(x -> occursin("post",x), readdir(n,join=true))
		out = filter(x -> !occursin(".1.",x),out)
		out = filter(x -> !occursin(".svg",x),out)
		iterMap = vcat(openAbcDistributions.(out)...);
		insertcols!(iterMap,1,:analysis => fill(analysis,size(iterMap,1)))
		insertcols!(iterMap,1,:b => fill(b,size(iterMap,1)))
		insertcols!(iterMap,1,:w => fill(aw,size(iterMap,1)))
		insertcols!(iterMap,1,:method => fill("ABC",size(iterMap,1)))

		dFiles = @. fName * "/div" * string(1:replicas) * ".tsv"
		function openTrues(d)
			tmp = Array(CSV.read(d,DataFrame))
			return(round.([tmp[3]/tmp[1],tmp[4]/tmp[1],(tmp[3]+tmp[4])/tmp[1]],digits=5))
		end

		trueAlphas = DataFrame(hcat(openTrues.(dFiles)...)',[:aw,:as,:a])

		insertcols!(trueAlphas,1,:analysis => fill(analysis,size(trueAlphas,1)))
		insertcols!(trueAlphas,1,:b => fill(b,size(trueAlphas,1)))
		insertcols!(trueAlphas,1,:w => fill(aw,size(trueAlphas,1)))
		insertcols!(trueAlphas,1,:method => fill("SLiM",size(trueAlphas,1)))

		push!(dmap,iterMap)
		push!(slimAlphas,trueAlphas)
	end

	dmap = vcat(dmap...)
	slimAlphas = vcat(slimAlphas...)
	
	dfAlphas  = vcat(dmap[:,1:end-2],slimAlphas)
	@rput dfAlphas

	p = R"""if(length(unique(dfAlphas[['w']])) > 1){
		dfAlphas[['w']] = factor(dfAlphas[['w']], labels = c('alpha[w]:0.1','alpha[w]:0.2','alpha[w]:0.3'))}else{dfAlphas[['w']] = factor(dfAlphas[['w']], labels = c('alpha[w]:0.1'))};
		if(length(unique(dfAlphas[['w']])) > 1){
		dfAlphas[['b']] = factor(dfAlphas[['b']], labels = c('B:0.2','B:0.4','B:0.8','B:0.999'))}else{dfAlphas[['w']] = factor(dfAlphas[['b']], labels = c('B:0.999'))}
		df = dfAlphas %>% arrange(b,w)
		df[["analysis"]] = factor(df[["analysis"]],levels=unique(df[["analysis"]]))
		d = reshape2::melt(df,id.vars=c("analysis","b","w","method"))
		p   = ggplot(d) + geom_boxplot(aes(x=variable,y=value,fill=method),alpha=0.8) + facet_wrap(~analysis,labeller=label_parsed,ncol=facetCols) + scale_fill_manual(values= c('#ab2710', '#e2bd9a', '#30504f'),name="Method") + colourSanMiguel() + theme_bw() + ggtitle($title) + xlab("Parameter") + ylab("Value");
		ggsave(p,filename=paste0($output,'.svg'),dpi=600)
		ggsave(p,filename=paste0($output,'.jpg'),dpi=600)"""
	return(dmap,dfAlphas)
end

