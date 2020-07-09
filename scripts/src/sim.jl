using Analytical, CSV,DataFrames, ProgressMeter

function summStats(param::Analytical.parameters,iter::Int64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)
	# @threads
	@showprogress for i in 1:iter
	# for i in 1:iter


		fac       = rand(-2:0.05:2)
		afac      = 0.184*(2^fac)
		bfac      = 0.000402*(2^fac)
		
		alTot     = rand(collect(0.05:0.01:0.4))
		lfac      = rand(collect(0.1:0.1:1))
		alLow     = round(alTot * lfac,digits=2)

		# println((thread=Threads.threadid(), iteration=i))
		bgsIter(param,afac,bfac,alTot,alLow,div,sfs,output,b,c)
	end
end

function bgsIter(param::Analytical.parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)

	for j in param.bRange
        # j = 0.999
        param.al = afac; param.be = bfac; 
        param.alLow = alLow; param.alTot = alTot; param.B = j

        Analytical.set_theta_f(param)
        theta_f = param.theta_f
        param.B = 0.999
        Analytical.set_theta_f(param)
        Analytical.setPpos(param)
        param.theta_f = theta_f
        param.B = j
        # x,y = Analytical.analyticalAlpha(param=param)
        x,y,z = Analytical.alphaByFrequencies(param,div,sfs,b,c)


        # println(z)

        Analytical.summaryStatistics(output, z)

		end
end


sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/noDemog/noDemog_0.4_0.2_0.999/sfs.tsv")))
sfs = sfs[:,2:end]
sfs = convert.(Int64,Analytical.cumulativeSfs(sfs))

sfsPos = [sfs[:,1] + sfs[:,2] sfs[:,1] + sfs[:,2]]
sfsNopos = [sfs[:,4] + sfs[:,2] sfs[:,4] + sfs[:,2]]

div = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/noDemog/noDemog_0.4_0.2_0.999/div.tsv")))
anDiv = [convert(Int64,sum(div)),convert(Int64,sum(div))]

# Set up model
adap = Analytical.parameters(N=500,n=661, gam_neg=-457, gL=10,gH=500,B=0.999,alTot=0.4,alLow=0.4)
#adap.nn=101
Analytical.binomOp(adap)

summStats(adap,1,anDiv,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/test",100,0.999)
summStats(adap,parse(Int,ARGS[1]),anDiv,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/" * ARGS[2] , 100,0.999)

