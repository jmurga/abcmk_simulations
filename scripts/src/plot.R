library(ggplot2)
library(data.table)

analysis = c("noDemog_0.4_0.1_0.2","noDemog_0.4_0.2_0.2","noDemog_0.4_0.3_0.2","noDemog_0.4_0.1_0.4","noDemog_0.4_0.2_0.4","noDemog_0.4_0.3_0.4","noDemog_0.4_0.1_0.8","noDemog_0.4_0.2_0.8","noDemog_0.4_0.3_0.8","noDemog_0.4_0.1_0.999","noDemog_0.4_0.2_0.999","noDemog_0.4_0.3_0.999")

analysis = c("tennesen_0.4_0.3_0.2","tennesen_0.4_0.3_0.4","tennesen_0.4_0.2_0.8","tennesen_0.4_0.3_0.8","tennesen_0.4_0.1_0.999","tennesen_0.4_0.2_0.999","tennesen_0.4_0.3_0.999")

plotPosteriorDensitiy = function(analysis,model,ne=500,path="/home/jmurga/mkt/202004/rawData/"){

	# analysis = c("noDemog_0.4_0.1_0.8","noDemog_0.4_0.2_0.8","noDemog_0.4_0.1_0.999","noDemog_0.4_0.2_0.999","noDemog_0.4_0.3_0.999")
	# model="noDemog";ne=500;path="/home/jmurga/mkt/202004/rawData/"
	out = list()
	density = list()
	segment = list()

	for(n in analysis){

		alphaWeak = unlist(strsplit(n,"_"))[3]
		bgs = unlist(strsplit(n,"_"))[4]

		df = fread(paste0(path,"/summStat/",model,"/",n,"/",n,"_posterior_1.tsv.gz"))
		df$alphaWeak= alphaWeak
		df$B = bgs
		df$analysis = n
		names(df) = c("Posterior alpha[w]","Posterior alpha[s]","Posterior alpha","alphaWeak","B","analysis")
		tmp = suppressWarnings(melt(df))
		
		density[[n]] = tmp
		div = fread(paste0(path,"/simulations/",model,"/ne",ne,"/",n,"/div.tsv"))

		aw = div$dw/div$di
		as = div$ds/div$di
		a = aw + as
		
		tmpSegment = data.table(x=c(aw,as,a),xend=c(aw,as,a),y=rep(0,3),yend=rep(Inf,3),c=c("#30504f", "#e2bd9a", "#ab2710"),g=c("True alpha[w]","True alpha[s]","True alpha"))
		tmpSegment$alphaWeak = alphaWeak
		tmpSegment$B = bgs
		tmpSegment$analysis = n
		segment[[n]] = tmpSegment

		p = ggplot(tmp) + 
			geom_density(aes(x=value,fill=variable),alpha=0.75) + 
			scale_fill_manual(values = paletteSanMiguel,labels=c(expression(paste("Posterior ",alpha[w])), expression(paste("Posterior ",alpha[s])),expression(paste("Posterior ",alpha)))) +
			geom_segment(data=tmpSegment,aes(x=x,xend=xend,y=y,yend=yend,colour=g),size=1, linetype=5) +
			scale_color_manual("True values",values=c("True alpha[w]"=paletteSanMiguel[1],"True alpha[s]"=paletteSanMiguel[2],"True alpha"=paletteSanMiguel[3]),labels=c(expression(paste("True ",alpha[w])), expression(paste("True ",alpha[s])),expression(paste("True ",alpha)))) + theme_bw() + labs(fill = "Posterior distributions",title=n,y="Posterior",x=expression(alpha))
		# ggsave(p,filename=paste0("/home/jmurga/mkt/202004/results/abc/",model,"/",n,".svg"),dpi=300)
	}

	dataPlot = rbindlist(density)	
	segment = rbindlist(segment)

	dataPlot$alphaWeak = factor(dataPlot$alphaWeak, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))
	segment$alphaWeak = factor(segment$alphaWeak, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))
	dataPlot$B = factor(dataPlot$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))
	segment$B = factor(segment$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))

	# segment$B = factor(segment$B, labels = c("B:0.8","B:0.999"))
	# dataPlot$B = factor(dataPlot$B, labels = c("B:0.999"))

	p2 = ggplot(dataPlot) + 
		geom_density(aes(x=value,fill=variable),alpha=0.75) + 
		scale_fill_manual(values = paletteSanMiguel,labels=c(expression(paste("Posterior ",alpha[w])), expression(paste("Posterior ",alpha[s])),expression(paste("Posterior ",alpha)))) +
		geom_segment(data=segment,aes(x=x,xend=xend,y=y,yend=yend,colour=g),size=1, linetype=5) +
		scale_color_manual("True values",values=c("True alpha[w]"=paletteSanMiguel[1],"True alpha[s]"=paletteSanMiguel[2],"True alpha"=paletteSanMiguel[3]),labels=c(expression(paste("True ",alpha[w])), expression(paste("True ",alpha[s])),expression(paste("True ",alpha)))) + 
		facet_grid(~analysis) + theme_bw() + facet_grid(B~alphaWeak,labeller=label_parsed) + labs(fill = "Posterior distributions")

	# Tables
	dataPlot$value = round(dataPlot$value,3)
	segment$x = round(segment$x,3)
	d1 = dataPlot %>% group_by(analysis,variable,alphaWeak,B) %>% summarize(q=paste0(round(mean(value),3)," [",quantile(value,c(0.1)),"-",quantile(value,0.9),"]"))
	d2 = reshape2::dcast(d1,analysis~variable)
	d3 = dcast(segment,analysis~g,value.var="x")
	d4 = merge(d2,d3,by='analysis') 
	d4 = d4[,c("analysis","True alpha[w]","Posterior alpha[w]","True alpha[s]","Posterior alpha[s]","True alpha","Posterior alpha")]

	tmp = do.call(rbind,strsplit(d4$analysis,"_"))[,3:4]
	d4$alphaWeak = tmp[,1]
	d4$B = tmp[,2]
	d4 = d4 %>% arrange(B,alphaWeak)

	# ggsave(p2,filename=paste0("/home/jmurga/mkt/202004/results/abc/",model,"/noDemog.svg"),dpi=300)
	out[['plot']] = p2
	out[['table']] = d4
	return(out)
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

# dfSegment = data.table(x=c(aw,as,a),xend=c(aw,as,a),y=rep(0,3),yend=rep(Inf,3))

# ggplot(d) + geom_density(aes(x=value,fill=variable,colour=variable),alpha=0.75) + fillSanMiguel() + colourSanMiguel() + geom_segment(data=dfSegment,aes(x=x,xend=xend,y=y,yend=yend),size=1,colour=dfSegment$c,linetype=2)            

r = plotPosterior(analysis=analysis,model="noDemog",pattern=NULL,nsample=10,software='R',ne=500)
ggsave(r$boxplot,filename="/home/jmurga/abcplots/noDemogABCRa.jpg",dpi=300)

r = plotPosterior(analysis=analysis,model="noDemog",pattern="postR",nsample=10,software='ABCtoolbox',ne=500)
ggsave(r$boxplot,filename="/home/jmurga/abcplots/noDemogABCRb.jpg",dpi=300)
ggsave(r$density,filename="/home/jmurga/abcplots/noDemogABCRb_density.jpg",dpi=300)

r = plotPosterior(analysis=analysis,model="noDemog",pattern="Best",nsample=10,software='ABCtoolbox',ne=500)
ggsave(r$boxplot,filename="/home/jmurga/abcplots/noDemogABCtoolbox.jpg",dpi=300)
ggsave(r$density,filename="/home/jmurga/abcplots/noDemogABCtoolbox_density.jpg",dpi=300)

#############################
r = plotPosterior(analysis=analysis,model="rescaled",pattern=NULL,nsample=50,software='R')
ggsave(r$boxplot,filename="/home/jmurga/abcplots/rescaledABCRa.svg")

r = plotPosterior(analysis=analysis,model="rescaled",pattern="postR",nsample=50,software='ABCtoolbox')
ggsave(r$boxplot,filename="/home/jmurga/abcplots/rescaledABCRb.svg")
ggsave(r$density,filename="/home/jmurga/abcplots/rescaledABCRb_density.svg")

r = plotPosterior(analysis=analysis,model="rescaled",pattern="Best",nsample=50,software='ABCtoolbox')
ggsave(r$boxplot,filename="/home/jmurga/abcplots/rescaledABCtoolbox.svg")
ggsave(r$density,filename="/home/jmurga/abcplots/rescaledABCtoolbox_density.svg")
#############################


analysis = c("tennesen_0.4_0.3_0.2","tennesen_0.4_0.3_0.4","tennesen_0.4_0.1_0.8","tennesen_0.4_0.2_0.8","tennesen_0.4_0.3_0.8","tennesen_0.4_0.1_0.999","tennesen_0.4_0.2_0.999","tennesen_0.4_0.3_0.999")

r = plotPosterior(analysis=analysis,model="tennesen",pattern=NULL,nsample=50,software='R')
ggsave(r$boxplot,filename="/home/jmurga/abcplots/tennesenABCRa.jpg",dpi=300)

r = plotPosterior(analysis=analysis,model="tennesen",pattern="postR",nsample=50,software='ABCtoolbox')
ggsave(r$boxplot,filename="/home/jmurga/abcplots/tennesenABCRb.jpg",dpi=300)
ggsave(r$density,filename="/home/jmurga/abcplots/tennesenABCRb_density.jpg",dpi=300)

r = plotPosterior(analysis=analysis,model="tennesen",pattern="Best",nsample=50,software='ABCtoolbox')
ggsave(r$boxplot,filename="/home/jmurga/abcplots/tennesenABCtoolbox.jpg",dpi=300)
ggsave(r$density,filename="/home/jmurga/abcplots/tennesenABCtoolbox_density.jpg",dpi=300)
