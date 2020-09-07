using Analytical, CSV,DataFrames, ProgressMeter

function summStats(param::Analytical.parameters,iter::Int64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)
    # @threads
    
        @showprogress for i in 1:iter
        # for i in 1:iter


                fac       = rand(-2:0.05:2)
                afac      = 0.184*(2^fac)
                bfac      = 0.000402*(2^fac)
                
                alTot     = rand(collect(0.01:0.01:0.4))
                # alTot     = rand(collect(0.05:0.05:0.4))
                # alLow     = round(rand(collect((alTot/10):(alTot/10):alTot)),digits=5)
                lfac      = rand(collect(0.1:0.1:0.9))
                alLow     = round(alTot * lfac,digits=5)
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

sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/"* ARGS[1] * "/" * ARGS[2] * "/sfs.tsv")))
# sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/tennesen/tennesen_0.4_0.2_0.999/sfs.tsv")))
sfs = sfs[:,2:end]
# sfs = convert.(Int64,Analytical.cumulativeSfs(sfs))

sfsPos   = sfs[:,1] + sfs[:,2]
sfsNopos = sfs[:,4] + sfs[:,2] 

divergence = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/"* ARGS[1] * "/" * ARGS[2] * "/div.tsv")))
# divergence = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/tennesen/tennesen_0.4_0.2_0.999/div.tsv")))
anDiv = [convert(Int64,sum(divergence[1:2]))]

# Set up modeld        = DataFrame(convert.(Int64,divergence))
pn       = sfs[1,1]
ps       = sfs[1,2]
rSfs     = Analytical.reduceSfs(sfs,100)'
alpha    = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'

# inputAbc = hcat(DataFrame(convert.(Int64,divergence[1:2]')),DataFrame([pn ps]),DataFrame(alpha),makeunique=true)
inputAbc = DataFrame(alpha)

CSV.write("/home/jmurga/mkt/202004/rawData/summStat/"* ARGS[1] * "/" * ARGS[2] * "/sfs"* ARGS[1] *".tsv", inputAbc, delim='\t',header=false);

adap = Analytical.parameters(N=731,n=661, gam_neg=-457, gL=10,gH=500,B=0.999,alTot=0.4,alLow=0.2)
#adap.nn=101
Analytical.binomOp!(adap)
summStats(adap,1,anDiv,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/test",100,0.9)

summStats(adap,parse(Int,ARGS[3]),anDiv,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/"* ARGS[1] * "/" * ARGS[2] * "/" * ARGS[2] * "_" * ARGS[4] , 100,0.9)
