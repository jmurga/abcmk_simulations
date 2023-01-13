using Distributed, JLD2, ProgressMeter, StatsBase
addprocs(8)
PATH = "/home/jmurga/mkt/202004/"

@everywhere using Analytical,DataFrames, CSV
Analytical.sourcePlotMapR(script="/home/jmurga/script.jl")
species = ["drosophila_melanogaster" "drosophila_santomea" "drosophila_simulans" "drosophila_teissieri" "drosophila_yakuba" "gorilla_gorilla" "homo_sapiens" "macaca_mulatta" "pan_troglodytes" "pongo_abelii"]

PATH="/home/jmurga/mkt/202004/rawData/nonModelOrganism"

out = []
for s in species
	analysis = PATH * "/" * s
	n = Int(countlines(analysis*"/sfs_"*s*".tsv")/2)
	adap = Analytical.parameters(N=1000,n=n)
	println(s,": ",n)

	h5file   = jldopen("/home/jmurga/mkt/202004/rawData/nonModelOrganism/rates_grapes.jld2")
	adap.dac = h5file["1000/"*string(n)*"/dac"]

	@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder=analysis,summstatSize=10^5,replicas=100,bootstrap=true);

	Analytical.ABCreg(analysis=analysis,replicas=100,P=5,S=size(adap.dac,1),tol=0.01,workers=20,abcreg="/home/jmurga/ABCreg/src/reg",parallel=true);

	df = Analytical.plotMap(analysis=analysis,weak=true,output = analysis * "/" * s * "_map.svg");
	insertcols!(df,1,:species => s)
	push!(out,df)
end

df = vcat(out...);
d = combine(groupby(df,:species),[:aw,:as,:a,:gamNeg,:shape] .=> mean)
d[:,2:end] = round.(d[:,2:end],digits=3);
d

CSV.write("abcmkPrimatesFlies.txt",d,delim='\t')

#=c("alpha","GammaExpo:negGmean","GammaExpo:negGshape","GammaExpo:posGmean","GammaExpo:pos_prop","prop_subst_sdel","prop_subst_wdel" ,"prop_subst_wadv","prop_subst_sadv")=#
#############################################

using Distributed, JLD2, ProgressMeter

PATH = "/home/jmurga/mkt/202004/"

@everywhere using Analytical,DataFrames, CSV
Analytical.sourcePlotMapR(script="/home/jmurga/script.jl")
adap = Analytical.parameters(N=1000,n=20)

h5file   = jldopen("/home/jmurga/mkt/202004/rawData/summStat/tbooker/tbooker.jld2")
adap.dac = h5file["1000/20/dac"]

@time summstat = Analytical.summaryStatsFromRates(param=adap,rates=h5file,analysisFolder="/home/jmurga/mkt/202004/rawData/summStat/tbooker/Nes10_pA0.0001/",summstatSize=10^5,replicas=100,bootstrap=true);

Analytical.ABCreg(analysis="/home/jmurga/mkt/202004/rawData/summStat/tbooker/Nes10_pA0.0001/",replicas=100,P=5,S=size(adap.dac,1),tol=0.01,workers=20,abcreg="/home/jmurga/ABCreg/src/reg",parallel=true);


df = Analytical.plotMap(analysis="/home/jmurga/mkt/202004/rawData/summStat/tbooker/Nes10_pA0.0001/",weak=false,output = "/home/jmurga/mkt/202004/rawData/summStat/tbooker/Nes10_pA0.0001/Nes10_pA0.0001_map.svg");
describe(df)

out = list()
trues = list()
for(i in list.dirs(recursive=F)){
	n = unlist(strsplit(i,"/"))[2]
	d = fread(paste0("/home/jmurga/mkt/202004/rawData/simulations/tbooker/",n,"/fix.tsv"))
	trues[[i]] = data.table(analysis=n,alpha=(d$dw/d$di))

    map = list()
    for (o in list.files(i,pattern="out",full.names=T)) {
        map[[o]] = apply(fread(o)[,2:5], 2, getmap) %>% t %>% as.data.table
    }
    tmp = rbindlist(map)
    tmp$analysis=n
    out[[n]] = tmp
}
out = rbindlist(out)
trues = rbindlist(trues)
d = melt(out[,c(2,5)],id.vars='analysis')0
#############################################

using Distributed, JLD2, ProgressMeter
addprocs(4)
PATH = "/home/jmurga/mkt/202004/"

@everywhere using Analytical,DataFrames, CSV,JLD2
@everywhere PATH = "/home/jmurga/mkt/202004/"
@everywhere include(PATH * "scripts/src/summaryParser.jl")


Analytical.sourcePlotMapR(script="/home/jmurga/script.jl")
adap = Analytical.parameters(N=1000,n=50)


analysis = ["twoepochs_1.15","twoepochs_1.25","twoepochs_1.5","twoepochs_5"]
n = size(analysis,1);

h5file = jldopen("/home/jmurga/twoepochs.jld2","r")

map(abcSummaries,analysis,fill(1000,n),fill(50,n),fill(h5file,n),fill(10^5,n),fill(100,n),fill(PATH * "rawData/summStat/twoepochs_bneck/",n),fill(false,n),fill("/home/jmurga/mkt/202004/rawData/simulations/twoepochs_bneck",n));

inputFolder = @. PATH * "rawData/summStat/twoepochs_bneck/" * analysis;
map(abcInference, inputFolder,fill(100,n),fill(5,n),fill(size(h5file["1000/50/dac"],1),n),fill(0.01,n),fill(20,n),fill(true,n));
df,dfAlphas = plotMapSimulations(inputFolder,100,inputFolder * "/" * "map.svg");