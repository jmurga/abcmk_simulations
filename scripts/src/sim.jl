using Analytical, CSV,DataFrames, ProgressMeter

function summStats(param::Analytical.parameters,iter::Int64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)
    # @threads
    
        @showprogress for i in 1:iter
        # for i in 1:iter

                fac       = rand(-2:0.05:2,2)
                afac      = 0.184*(2^fac[1])
                bfac      = 0.000402*(2^fac[2])
                
                alTot     = rand(collect(0.01:0.01:0.6))
                lfac      = rand(collect(0.1:0.1:0.9))
                alLow     = round(alTot * lfac,digits=5)

                bgsIter(param,afac,bfac,alTot,alLow,div,sfs,output,b,c)
        end
end

function bgsIter(param::Analytical.parameters,afac::Float64,bfac::Float64,alTot::Float64,alLow::Float64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)

        for j in param.bRange
                # j = 0.999
                param.al = afac; param.be = bfac; 
                param.alLow = alLow; param.alTot = alTot; param.B = j

                Analytical.set_theta_f!(param)
                theta_f = param.theta_f
                param.B = 0.999
                Analytical.set_theta_f!(param)
                Analytical.setPpos!(param)
                param.theta_f = theta_f
                param.B = j
                #x,y,z = Analytical.alphaByFrequencies(param,d,sum(sfs[:,1:2],dims=2),100,0.999)

                # x,y = Analytical.analyticalAlpha(param=param)
                x,y,z = Analytical.alphaByFrequencies(param,div,sfs,b,c)
                # println(z)

                Analytical.summaryStatistics(output, DataFrame(z))

        end
end

sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/"* ARGS[1] * "/" * ARGS[2] * "/sfs.tsv")))
sfs = sfs[:,2:end]

sCumu = convert.(Int64,Analytical.cumulativeSfs(sfs))

sfsPos   = sfs[:,1] + sfs[:,2]
sfsNopos = sfs[:,4] + sfs[:,2] 

divergence = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/"* ARGS[1] * "/" * ARGS[2] * "/div.tsv")))
d = [convert(Int64,sum(divergence[1:2]))]

# rSfs     = Analytical.reduceSfs(sfs,100)'
# alpha    = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'
rSfsCumu    = Analytical.cumulativeSfs(Analytical.reduceSfs(sfs,100)')
rSfs   = Analytical.reduceSfs(sfs,100)'
# rSfsCumu    = Analytical.cumulativeSfs(sfs)
alphaCumu   = @. round(1 - divergence[2]/divergence[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)'
alphaReduced   = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'

# inputAbc = hcat(DataFrame(convert.(Int64,divergence[1:2]')),DataFrame([pn ps]),DataFrame(alpha),makeunique=true)
inputAbc = DataFrame(alphaReduced)

CSV.write("/home/jmurga/mkt/202004/rawData/summStat/"* ARGS[1] * "/" * ARGS[2] * "/sfs"* ARGS[1] *".tsv", inputAbc, delim='\t',header=false);

adap = Analytical.parameters(N=parse(Int,ARGS[3]),n=parse(Int,ARGS[4]), gam_neg=-457, gL=10,gH=500,B=0.999,alTot=0.4,alLow=0.2)
#adap.nn=101
#adap.nn = adap.nn + 1
Analytical.binomOp!(adap)
summStats(adap,1,d,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/test",100,0.999)

summStats(adap,parse(Int,ARGS[5]),d,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/" * ARGS[1] * "/" * ARGS[2] * "/" * ARGS[2] * "_" * ARGS[6] , 100,0.999)
