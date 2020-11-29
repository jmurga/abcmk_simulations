library(ggplot2)
library(data.table)

analysis = c("noDemog_0.4_0.1_0.2","noDemog_0.4_0.2_0.2","noDemog_0.4_0.3_0.2","noDemog_0.4_0.1_0.4","noDemog_0.4_0.2_0.4","noDemog_0.4_0.3_0.4","noDemog_0.4_0.1_0.8","noDemog_0.4_0.2_0.8","noDemog_0.4_0.3_0.8","noDemog_0.4_0.1_0.999","noDemog_0.4_0.2_0.999","noDemog_0.4_0.3_0.999")

plotPosterior = function(analysis,ne){

	alphas = list()
	out = list()
	density = list()
	output = list()

	model = unlist(strsplit(analysis,"_"))[1]

	for(n in analysis){
		print(n)		
		bgs = unlist(strsplit(n,"_"))[4]
		aw = unlist(strsplit(n,"_"))[3]

		sim= paste0("/home/jmurga/mkt/202004/rawData/simulations/",model,"/ne",ne,"/",n)
		ss = paste0("/home/jmurga/mkt/202004/rawData/summStat/",model,"/",n)

		l <- vector("list", length = length(list.files(ss,pattern='posterior')))
		a <- list()
		tmp = fread(paste0(sim,"/divTable.tsv"))
		# tmp = fread(paste0(sim,"/alphas.tsv"))
		# tmp$analysis = n
		# tmp$method = 'simulations'
		# tmp$B = bgs
		# tmp$alphaWeak = aw
		# names(tmp) = c('alphaW','alphaS','alpha','analysis',"method","B","alphaWeak")	

		for(i in list.files(ss,pattern='posterior')){
			df <- fread(paste0(ss,"/",i))
			df$analysis = n
			density[[n]][[i]] = df
			l[[i]] <-  transpose(data.table(colMeans(x=df[,1:3], na.rm = TRUE)))

			alpha = tmp[sample(.N,500,replace=F)] %>% summarize_all(sum) %>% mutate(alphaW=dw/di,alphaS=ds/di,alpha=(ds+dw)/di) %>% select(alphaW,alphaS,alpha)

			a[[i]] = alpha
			a[[i]]$analysis = n
			a[[i]]$method = "Simulations"
			a[[i]]$B = bgs
			a[[i]]$alphaWeak = aw
			names(a[[i]]) = c('alphaW','alphaS','alpha','analysis',"method","B","alphaWeak")
		}
		
		alphas[[n]] = rbindlist(a)
		# alphas[[n]] = tmp
		density[[n]] = rbindlist(density[[n]])
		density[[n]]$B = bgs
		density[[n]]$alphaWeak = aw
		names(density[[n]]) <- c('alphaW','alphaS','alpha','analysis','B','alphaWeak')

		tmp = rbindlist(l)
		tmp$analysis = n
		tmp$method = "ABC"
		tmp$B = bgs
		tmp$alphaWeak = aw
		names(tmp) <- c('alphaW','alphaS','alpha','analysis',"method","B","alphaWeak")

		out[[n]] = tmp
	}

	alphas = rbindlist(alphas)
	out = rbindlist(out)

	dataPlot = rbind(alphas,out)

	dataMelt = melt(dataPlot,id.vars=c("analysis","method","B","alphaWeak"))
	dataMelt$analysis = factor(dataMelt$analysis,levels=analysis)


	# dataPlot %>% group_by(method,analysis) %>% summarise(	output = list()alphaW = quantile(alphaW, c(0.25,0.75)), q = c(0.25, 0.75),alphaS = quantile(alphaS, c(0.25,0.75)), q = c(0.25, 0.75))

	# out = c(dataPlot,p)
	# return(dataPlot)

	# ggsave(p,file="/home/jmurga/noDemog_ABC.svg")

	d = melt(rbindlist(density))
	d$analysis = factor(d$analysis,levels=analysis)
	# d$analysis = factor(d$analysis,levels=analysis)

	output[['out']] = dataMelt
	output[['density']] = d

	return(output)
}

plotSimulatedAlphas = function(f,title,output="/home/jmurga/mkt/202004/results/simulations/alphasSimulations/"){

	df = fread(f)
	df$alphaW = factor(df$alphaW, labels = c(
    "alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))	
    df$B = factor(df$B, labels = c(
    "B:0.2","B:0.4","B:0.8","B:0.999"))

	p = ggplot(df, aes(x=f,y=alphas,color=sfs)) + 
		geom_line() + 
		facet_grid(B~alphaW,labeller=label_parsed) + 
		ylim(-0.5,0.5) + 
		scale_color_manual(values=c('black','#ab2710'),name = "SFS", labels = c("All alleles","Neutral + deleterious")) + 
		theme_bw() + 
		ylab(expression(alpha)) + 
		xlab('Derived Allele Count') + 
		theme(legend.position="bottom",legend.text=element_text(size=14),plot.title=element_text(hjust=0.5,face='bold')) +
		ggtitle(title);p

		ggsave(p,filename=paste0("/home/jmurga/mkt/202004/results/simulations/alphasSimulations/",output,".svg"))
		ggsave(p,filename=paste0("/home/jmurga/mkt/202004/results/simulations/alphasSimulations/",output,".jpg"),dpi=600)

}

r = plotPosterior(analysis)

df=r$out

df$alphaWeak = factor(df$alphaWeak, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))	
df$B = factor(df$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))
p = ggplot(df,aes(x=variable,y=value,fill=method)) + geom_boxplot() + scale_fill_manual(values=c("#ab2710","#e2bd9a")) + facet_grid(~analysis)+ theme_bw() + facet_grid(B~alphaWeak,labeller=label_parsed) +scale_x_discrete(labels = c(expression(alpha[w]),expression(alpha[s]),expression(alpha)))

ggsave(p,filename="/home/jmurga/tennesenABC.svg")

ggsave(p,filename="/home/jmurga/mkt/202004/results/abc/boxplots/noDemogABC.svg")
ggsave(p,filename="/home/jmurga/mkt/202004/results/abc/boxplots/noDemogABC.jpg",dpi=600)
fwrite(df,"/home/jmurga/mkt/202004/results/abc/boxplots/noDemogABC.tsv")

fwrite(r$density,"/home/jmurga/mkt/202004/results/abc/boxplots/noDemogDensity.tsv")


r$density$alphaWeak = factor(r$density$alphaWeak, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))	
r$density$B = factor(r$density$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))
pD <- ggplot(r$density,aes(x=value,fill=variable)) + geom_density() + scale_fill_manual(values=c("#30504f","#e2bd9a","#ab2710")) + facet_grid(~analysis)+ theme_bw() + facet_grid(B~alphaWeak,labeller=label_parsed) +scale_x_discrete(labels = c(expression(alpha[w]),expression(alpha[s]),expression(alpha)));pD

# ggsave(pD,filename="/home/jmurga/noDemog_density.svg")
