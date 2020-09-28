using Analytical, CSV,DataFrames, ProgressMeter

function summStats(param::Analytical.parameters,iter::Int64,div::Array,sfs::Array,output::String,b::Int64,c::Float64)
    # @threads

        @showprogress for i in 1:iter
        # for i in 1:iter

            fac       = rand(-2:0.05:2,2)
            afac      = 0.184*(2^fac[1])
            bfac      = 0.000402*(2^fac[2])

            alTot     = rand(collect(0.05:0.01:0.4))
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


model, path, popSize, nSize, bins, iterations, iterFile = ARGS
#model, path, bins= "noDemog", "noDemog_0.4_0.2_0.999","100"

sfs = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/"* model * "/" * path * "/sfs.tsv")))
sfs = sfs[:,2:end]

sCumu = convert.(Int64,Analytical.cumulativeSfs(sfs))

sfsPos   = sfs[:,1] + sfs[:,2]
sfsNopos = sfs[:,4] + sfs[:,2]

divergence = convert(Array,DataFrame!(CSV.File("/home/jmurga/mkt/202004/rawData/simulations/"* model * "/" * path * "/div.tsv")))
d = [convert(Int64,sum(divergence[1:2]))]

# rSfs     = Analytical.reduceSfs(sfs,100)'
# alpha    = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'
rSfsCumu    = Analytical.cumulativeSfs(Analytical.reduceSfs(sfs,parse(Int,bins))')
rSfs   = Analytical.reduceSfs(sfs,parse(Int,bins))'
# rSfsCumu    = Analytical.cumulativeSfs(sfs)
alphaCumu   = @. round(1 - divergence[2]/divergence[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)'
alphaReduced   = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'

# inputAbc = hcat(DataFrame(convert.(Int64,divergence[1:2]')),DataFrame([pn ps]),DataFrame(alpha),makeunique=true)
#inputAbc1 = DataFrame(alphaReduced)
inputAbc2 = DataFrame(alphaCumu)

#CSV.write("/home/jmurga/mkt/202004/rawData/summStat/"* model * "/" * path * "/sfs"* model * "_" * bins * ".tsv", inputAbc1, delim='\t',header=false);
CSV.write("/home/jmurga/mkt/202004/rawData/summStat/"* model * "/" * path * "/sfs"* model * "_" * bins * ".tsv", inputAbc2, delim='\t',header=false);

adap = Analytical.parameters(N=parse(Int,popSize),n=parse(Int,nSize), gam_neg=-457, gL=10,gH=500,B=0.999,alTot=0.4,alLow=0.2,Lf=2*10^5)
#adap.nn=101
#adap.nn = adap.nn + 1
Analytical.binomOp!(adap)
summStats(adap,1,d,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/test",parse(Int,bins),0.999)

summStats(adap,parse(Int,iterations),d,sfsPos,"/home/jmurga/mkt/202004/rawData/summStat/" * model * "/" * path * "/" * path * "_" * iterFile , parse(Int,bins),0.999)
