using Analytical, CSV,DataFrames, ProgressMeter

model, path, popSize, nSize, bins, iterations, iterFile = ARGS

function readSfs(models,path="/home/jmurga/mkt/202004/rawData/simulations/")
    
    dirs = path .* split.(models,'_')[1] .* '/' .*models;

    s = []; d = []; a = [];

    for f in eachindex(dirs) 
        
        sfs = convert(Array,DataFrame!(CSV.File(dirs[f] * "/sfs.tsv")));   
        sCumu = convert.(Int64,Analytical.cumulativeSfs(sfs[:,2:end]))
        sfsPos   = sCumu[:,1] + sCumu[:,2]

        div1 = convert(Array,DataFrame!(CSV.File(dirs[f] * "/div.tsv")));
        div2 = convert(Int64,sum(div1[1:2]))

        rSfsCumu    = Analytical.cumulativeSfs(Analytical.reduceSfs(sCumu,parse(Int,bins))')
        # rSfsCumu    = Analytical.cumulativeSfs(sfs)
        alphaCumu   = @. round(1 - div1[2]/div1[1] * rSfsCumu[:,1]/rSfsCumu[:,2],digits=5)'

        push!(s, sfsPos)
        push!(d, divergence)
        push!(a, alphaCumu)
    end
    
    s = reduce(hcat,s)

    return(s,d,a)
end


# rSfs     = Analytical.reduceSfs(sfs,100)'
# alpha    = @. round(1 - divergence[2]/divergence[1] * rSfs[:,1]/rSfs[:,2],digits=5)'


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
