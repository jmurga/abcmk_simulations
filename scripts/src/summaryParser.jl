function readSfsDiv(analysis,path="/home/jmurga/mkt/202004/rawData/simulations/noDemog/")

	sfs = convert(Array,CSV.read(path * analysis * "/sfs.tsv"))
	div = convert.(Int64,convert(Array,CSV.read(path * analysis * "/div.tsv")))

	s = convert.(Int64,Analytical.cumulativeSfs(sfs[:,2:end])[:,1:2])

	return s, div

end

function makepriorsFiles(;df,alpha,samples,nfiles,analysis,bins,nthreads,path="/home/jmurga/mkt/202004/rawData/summStat/noDemog",parallel="/home/jmurga/.conda/envs/abcmk/bin/parallel",gzip="/usr/bin/gzip")

	dir = `mkdir -p $path/bins$bins/$analysis`;
	run(dir);

	# Writing sampled summaries
	for i in 1:nfiles
		idx = sample(eachrow(df).iter,samples;replace=false)
		tmp = df[idx,:]

		CSV.write(path * "/" * "bins" * string(bins) * "/" * analysis * "/" * analysis * "_" * string(i) * ".tsv", DataFrames.DataFrame(tmp), delim='\t', append=false,header=false)

	end

	# Gzip summaries
	f = path * "/bins" * string(bins) * "/" * analysis
	f = f .* "/" .* analysis .* "_" .* string.(collect(1:nfiles)) .* ".tsv"

	gz = `$parallel -j$nthreads $gzip -f '{''}' ::: $f`;
	run(gz);

	# Writting and gzip alpha
	CSV.write(path * "/" * "bins" * string(bins) * "/" * analysis * "/alpha_" * analysis * ".tsv", DataFrames.DataFrame(alpha), delim='\t', append=false,header=false)
	fAlpha = path * "/" * "bins" * string(bins) * "/" * analysis * "/alpha_" * analysis * ".tsv"
	gz = `$gzip $fAlpha`;
	run(gz);
end

function runSummary(;param,nsamples,iterations,nfiles,bins,nthreads,analysis)
	sfs, d = readSfsDiv(analysis);

	# Making alpha
	sfs, d = readSfsDiv(analysis);
	rSfs = Analytical.reduceSfs(sfs,bins)'
	α    = @. round(1 - (d[2]/d[1] * rSfs[:,1]/rSfs[:,2]),digits=5)'
	# Making summary statistics
	@time df = Analytical.summaryStats(param = adap,alpha = 0.4, divergence = [sum(d)], sfs = sum(sfs,dims=2),bins = bins ,iterations = iterations);
	makepriorsFiles(df=df,alpha = α,samples=nsamples, nfiles=nfiles,analysis=analysis,bins=bins,nthreads=nthreads)

	return df
end
