library(ggplot2)

analysis = c("noDemog_0.4_0.1_0.2","noDemog_0.4_0.2_0.2","noDemog_0.4_0.3_0.2","noDemog_0.4_0.1_0.4","noDemog_0.4_0.2_0.4","noDemog_0.4_0.3_0.4","noDemog_0.4_0.1_0.8","noDemog_0.4_0.2_0.8","noDemog_0.4_0.3_0.8","noDemog_0.4_0.2_0.999","noDemog_0.4_0.1_0.999","noDemog_0.4_0.3_0.999")

analysis = c("noDemog_0.4_0.1_0.2")
bins=50
alphas = list()
output = list()
density = list()
for(n in analysis){
	print(n)
	sim= paste0("/home/jmurga/mkt/202004/rawData/simulations/noDemog/",n)
	ss = paste0("/home/jmurga/mkt/202004/rawData/summStat/noDemog/bins",bins,"/",n)
	alphas[[n]] = fread(paste0(sim,"/alphas.tsv"))
	alphas[[n]]$analysis = n
	alphas[[n]]$method = "simulation"
	# alphas[[n]] = alphas[[n]][sample(.N,100,replace=F)]
    
    names(alphas[[n]]) = c('alphaW','alphaS','alpha','analysis',"method")
	
	l <- vector("list", length = length(list.files(ss,pattern='0.tangent')))
	for(i in list.files(ss,pattern='0.tangent')){
		df <- fread(paste0(ss,"/",i))
		df$analysis = n
		density[[n]][[i]] = df
		l[[i]] <-  transpose(data.table(colMeans(x=df[,1:3], na.rm = TRUE)))
	}

	density[[n]] = rbindlist(density[[n]])
	names(density[[n]]) <- c('alphaW','alphaS','alpha','analysis')

	tmp = rbindlist(l)
	tmp$analysis = n
	tmp$method = "abc"
	names(tmp) <- c('alphaW','alphaS','alpha','analysis',"method")
	
	output[[n]] = tmp
}

alphas = rbindlist(alphas)
output = rbindlist(output)

dataPlot = rbind(alphas,output)
dataMelt = melt(dataPlot)
# dataMelt$analysis = factor(dataMelt$analysis,levels=c("noDemog_0.4_0.1_0.999","noDemog_0.4_0.2_0.999","noDemog_0.4_0.3_0.999","noDemog_0.4_0.3_0.8"))
dataMelt$analysis = factor(dataMelt$analysis,levels=c("noDemog_0.4_0.2_0.2"))

# result = rbind(df,alphas[sample(.N,100,replace=F)]) %>% melt
p <- ggplot(dataMelt,aes(x=variable,y=value,fill=method)) + geom_boxplot() + scale_fill_manual(values=c("#ab2710","#e2bd9a")) + facet_wrap(~ analysis, ncol=3)

ggsave(p,file="/home/jmurga/noDemog0.2.png",dpi=300)
ggsave(p,file="/home/jmurga/noDemog.png",dpi=300)

d = melt(rbindlist(density))
# d$analysis = factor(d$analysis,levels=c("noDemog_0.4_0.1_0.999","noDemog_0.4_0.2_0.999","noDemog_0.4_0.3_0.999","noDemog_0.4_0.3_0.8"))
d$analysis = factor(d$analysis,levels=c("noDemog_0.4_0.2_0.2"))

pD <- ggplot(d,aes(x=value,fill=variable)) + geom_density() + scale_fill_manual(values=c("#30504f","#e2bd9a","#ab2710")) + facet_wrap(~ analysis,ncol=3)

dataAlpha = melt(alphas)
pA = ggplot(dataAlpha,aes(x=value,fill=variable)) + geom_density() + scale_fill_manual(values=c("#30504f","#e2bd9a","#ab2710")) + facet_wrap(~ analysis,ncol=3)

ggsave(pA,file="/home/jmurga/noDemogDensityReal.png",dpi=300)








# rm jjob.txt;for i in {1..15};do printf "julia /home/jmurga/mkt/202004/scripts/src/sim.jl noDemog_0.4_0.1_0.8 39300 ${i}\n" >> jjob.txt;done

# for i in {1..20};do cat noDemog_0.4_0.1_0.999_${i}_1.tsv >> slim.tsv;done

# rm jjob.txt;for i in {1..20};do printf "julia /home/jmurga/mkt/202004/scripts/src/sim.jl noDemog_0.4_0.3_0.999 3000 ${i}\n" >> jjob.txt;done
# for i in {1..20};do cat noDemog_0.4_0.3_0.999_${i}_1.tsv >> slim.tsv;done

# rm jjob.txt;for i in {1..20};do printf "julia /home/jmurga/mkt/202004/scripts/src/sim.jl noDemog_0.4_0.1_0.999 3000 ${i}\n" >> jjob.txt;done
# for i in {1..20};do cat noDemog_0.4_0.1_0.999_${i}_1.tsv >> slim.tsv;done


# rm jjob.txt;for i in {1..20};do printf "julia /home/jmurga/mkt/202004/scripts/src/sim.jl noDemog_0.4_0.2_0.8 3000 ${i}\n" >> jjob.txt;done
# for i in {1..20};do cat noDemog_0.4_0.2_0.8_${i}_1.tsv >> slim.tsv;done

# rm jjob.txt;for i in {1..20};do printf "julia /home/jmurga/mkt/202004/scripts/src/sim.jl noDemog_0.4_0.3_0.8 3000 ${i}\n" >> jjob.txt;done
# for i in {1..20};do cat noDemog_0.4_0.3_0.8_${i}_1.tsv >> slim.tsv;done

