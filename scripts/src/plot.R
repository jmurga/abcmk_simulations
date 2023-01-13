library(ggplot2)
library(data.table)

mapDistribution = function(case,path="/home/jmurga/mkt/202004/rawData/summStat/"){

	lmean = list()
	# lmode = list()
	lmap  = list()
	lnegative  = list()
	ltrues  = list()

	output = list()
	for(n in list.dirs(paste0(path,case),recursive=F)){

		analysis = unlist(strsplit(n,"/"))
		analysis = analysis[length(analysis)]
		print(analysis)

		b = as.numeric(unlist(strsplit(analysis,"_"))[4])
		w = as.numeric(unlist(strsplit(analysis,"_"))[3])
		a = unlist(strsplit(analysis,"_"))[1]
		
		out = list.files(n,pattern='out',full.names=T)
		# out = list.files(pattern='out',full.names=T)
		out = out[!grepl('.1.post.gz',out)]
		out = out[!grepl('tangent',out)]
		

		#########
		# s = fread(paste0("/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/",analysis,"/sfs.tsv"))
		# d = fread(paste0("/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/",analysis,"/div.tsv"))
		# sfs = cumulative(s)
		# a = alpha(sfs,d,pw=T)[700]
		# d = data.table(a*w/0.4,a * (1-w/0.4),a,analysis,w,b)
		#########	
		d = data.table(readDiv(analysis),analysis,w,b)
		
		# d$aw = d$a * w/0.4
		# d$as = d$a * (1-w/0.4)
		d = melt(d,id.vars=c('analysis','b','w'))
		d$group=c("#30504f", "#e2bd9a", "#ab2710")
		ltrues[[n]] = d

		for(i in out){
			df = fread(i)
			neg = df[,4:5]
			df = df[,1:3]
			names(df) = c('aw','as','a')
			names(neg) = c('gamNeg','shape')
			mn = data.table(t(apply(df,2,mean)),analysis,b,w)
			# md = data.table(t(apply(df,2,Mode)),analysis,b,w)
			map = data.table(t(apply(df,2,getmap)),analysis,b,w)
			mneg = data.table(t(apply(neg,2,getmap)),analysis,b,w)

			lmean[[i]] = mn
			# lmode[[i]] = md
			lmap[[i]] = map 
			lnegative[[i]] = mneg
		}
	}

	pnegative = rbindlist(lnegative)
	tableNeg = pnegative %>% group_by(analysis,b,w) %>% summarise(gamNeg = paste0(round(mean(gamNeg),3)," [",round(quantile(gamNeg,0.05),3),",",round(quantile(gamNeg,0.95),3),"]"),shape = paste0(round(mean(shape),3)," [",round(quantile(shape,0.05),3),",",round(quantile(shape,0.95),3),"]"))

	pnegative = pnegative %>% arrange(b,w)
	pnegative$analysis = factor(pnegative$analysis,levels=unique(pnegative$analysis))
	pnegative = melt(pnegative,id.vars=c('analysis','b','w'))
	pg = ggplot(pnegative[variable == 'gamNeg']) + geom_density(aes(value,fill=variable),alpha=0.5) + facet_wrap(~analysis,ncol=3) + fillSanMiguel() + geom_vline(xintercept=-457)

	ps =  ggplot(pnegative[variable == 'shape']) + geom_density(aes(value,fill=variable),alpha=0.5) + facet_wrap(~analysis,ncol=3) + fillSanMiguel() + geom_vline(xintercept=0.184)

	trues=rbindlist(ltrues)
	trues = trues %>% arrange(b,w)
	trues$analysis = factor(trues$analysis,levels=unique(trues$analysis))

	pmap = rbindlist(lmap)
	pmap = pmap %>% arrange(b,w)
	pmean = rbindlist(lmean)
	pmean = pmean %>% arrange(b,w)

	pmap$analysis= factor(pmap$analysis,levels=unique(pmap$analysis))
	pmean$analysis= factor(pmean$analysis,levels=unique(pmean$analysis))
	df = melt(pmap,id.vars=c('analysis','b','w'))
	df2 = melt(pmean,id.vars=c('analysis','b','w'))

	p = ggplot(df) + geom_density(aes(x=value,fill=variable),alpha=0.5) + facet_wrap(~analysis,ncol=3) + fillSanMiguel() + geom_segment(data=trues,aes(x=value,xend=value,y=0,yend=Inf,colour=variable),size=1, linetype=5) + colourSanMiguel()

	p2 = ggplot(df2) + geom_density(aes(x=value,fill=variable),alpha=0.5) + facet_wrap(~analysis,ncol=3) + fillSanMiguel() + geom_segment(data=trues,aes(x=value,xend=value,y=0,yend=Inf,colour=variable),size=1, linetype=5) + colourSanMiguel()

	diff1 = pmap %>% group_by(analysis,b,w) %>% summarise_all(mean) %>% as.data.table
	diff2 = dcast(trues,analysis+b+w~variable) %>% as.data.table
	diff  = cbind(diff1[,1:3],round(diff1[,4:6] - diff2[,4:6],5))

	pmap$case=case
	diff$case=case
	trues$case=case

	output[['diff']]      = diff
	output[['negative']]  = tableNeg
	output[['trues']]     = trues
	output[['map']]       = map
	output[['plot']]      = p
	output[['plot_mean']] = p2
	
	return(output)
}

mapDistributionSlim = function(case,path="/home/jmurga/mkt/202004/rawData/summStat/"){

	case="noDemog_bgs_v2"
	path="/home/jmurga/mkt/202004/rawData/summStat/"
	# n= list.dirs(paste0(path,case),recursive=F)[4]

	lmean      = list()
	lnegative  = list()
	lmap       = list()
	ltrues     = list()
	slimAlphas = list()
	output     = list()
	for(n in list.dirs(paste0(path,case),recursive=F)){

		analysis = unlist(strsplit(n,"/"))
		analysis = analysis[length(analysis)]
		print(analysis)

		b = as.numeric(unlist(strsplit(analysis,"_"))[4])
		w = as.numeric(unlist(strsplit(analysis,"_"))[3])
		a = unlist(strsplit(analysis,"_"))[1]
		
		out = list.files(n,pattern='log',full.names=T)
		# out = list.files(pattern='out',full.names=T)
		out = out[!grepl('tangent',out)]
		out = out[!grepl('.1.post.gz',out)]

		for(i in out){
			df = fread(i)
			neg = df[,4:5]
			df = df[,1:3]
			names(df) = c('aw','as','a')
			names(neg) = c('gamNeg','shape')
			mn = data.table(t(apply(df,2,mean)),analysis,b,w)
			# md = data.table(t(apply(df,2,Mode)),analysis,b,w)
			map = data.table(t(apply(df,2,getmap)),analysis,b,w)
			mneg = data.table(t(apply(neg,2,getmap)),analysis,b,w)

			lmean[[paste0(n,i)]] = mn
			lmap[[paste0(n,i)]] = map 
			lnegative[[paste0(n,i)]] = mneg
		}

		for(i in 1:100){
			print(i)
			sfs = fread(paste0("/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/",analysis,"/sfs",i,".tsv"),header=F)
			names(sfs) = c('f','pi','p0','pw')
			sfs$pi_nopos = sfs$pi - sfs$pw
			scumu = cumulative(sfs)

			div = fread(paste0("/home/jmurga/mkt/202004/rawData/simulations/noDemog/ne500/",analysis,"/div",i,".tsv"),header=F)
			tmp = (div[3]+div[4])/div[1]
			daf = scumu[,c(1,5,3)]
			names(daf) = c('daf','Pi','P0')
			names(div) = c('Di','D0','Dw','Ds')
			div$mi = 10^8*0.75
			div$m0 = 10^8*0.25

			# tmp = 1 - (div$D0/div$Di * scumu$pi_nopos/scumu$p0)[700]
			# tmp = 1 - (div$D0/div$Di * scumu$pi_nopos/scumu$p0)[700]
			tmp = aMKT(daf,div[,c(1,2,5,6)],xhigh=0.75)$alphaCorrected$alphaAsymptotic
			slimAlphas[[paste0(analysis,i)]] = data.table(aw= tmp * (w/0.4),as=tmp * (1-(w/0.4)),a = tmp,analysis=analysis,b=b,w=w)
		}
	}

	# Join data
    trues          = rbindlist(slimAlphas)
    trues          = trues %>% arrange(b,w)
    trues$analysis = factor(trues$analysis,levels=unique(trues$analysis))
    trues$method   = 'slim'

    pmap           = rbindlist(lmap)
    pmap           = pmap %>% arrange(b,w)
    pmap$method    = 'abc'
	
    pmean          = rbindlist(lmean)
    pmean          = pmean %>% arrange(b,w)
    pmean$method   = 'abc'
	
    pmap$analysis  = factor(pmap$analysis,levels=unique(pmap$analysis))
    pmean$analysis = factor(pmean$analysis,levels=unique(pmean$analysis))

    # To plot
    df  = rbind(pmap,trues)
    d   = melt(df,id.vars=c('analysis','b','w','method'))
    df2 = melt(pmean,id.vars=c('analysis','b','w'))

	# p = ggplot(df) + geom_density(aes(x=value,fill=variable),alpha=0.5) + facet_wrap(~analysis,ncol=3) + fillSanMiguel() + geom_segment(data=trues,aes(x=value,xend=value,y=0,yend=Inf,colour=variable),size=1, linetype=5) + colourSanMiguel()

	# p2 = ggplot(df2) + geom_density(aes(x=value,fill=variable),alpha=0.5) + facet_wrap(~analysis,ncol=3) + fillSanMiguel() + geom_segment(data=trues,aes(x=value,xend=value,y=0,yend=Inf,colour=variable),size=1, linetype=5) + colourSanMiguel()

	# diff1 = pmap %>% group_by(analysis,b,w) %>% summarise_all(mean) %>% as.data.table
	# diff2 = dcast(trues,analysis+b+w~variable) %>% as.data.table
	# diff  = cbind(diff1[,1:3],round(diff1[,4:6] - diff2[,4:6],5))

	# pmap$case=case
	# diff$case=case
	# trues$case=case

 #    output[['diff']]  = diff
 #    output[['trues']] = trues
 #    output[['map']]   = map
 #    output[['plot']]  = p
 #    output[['plot_mean']]  = p2
	return(output)
}

plotPosteriorDensity = function(path="/home/jmurga/mkt/202004/rawData/summStat/noDemog",global=FALSE){
	
	# model="noDemog";ne=500;path="/home/jmurga/mkt/202004/rawData/";global=F;fixed=F
	# model="isolation";ne=NULL;path="/home/jmurga/mkt/202004/rawData";global=F;fixed=F
	
	out = list()
	density = list()
	rejection = list()
	segment = list()
	modes = list()
	analysis = list.dirs(path,recursive=F)


	for(n in analysis){
		print(n)
		tmp   = unlist(strsplit(n,"_"))
		model = unlist(strsplit(tmp[1],"/"))
		model = model[length(model)]
		alphaWeak = tmp[3]
		bgs  = tmp[4]


		nd = unlist(strsplit(n,"/"))[length(unlist(strsplit(n,"/")))]

		if(model == 'isolation'){
			div = fread(paste0("/home/jmurga/mkt/202004/rawData/simulations/",model,"/",nd,"/div.tsv"))
			
			# df1 = fread(paste0(path,"/summStat/",model,"/",n,"/out.0.log.post.gz"))
			# df2 = fread(paste0(path,"/summStat/",model,"/",n,"/neuralnet_1.tsv.gz"))
			df1 = fread(paste0(path,"/",nd,"/","out_1.0.post.gz"))
			df2 = fread(paste0(path,"/",nd,"/","out_1.0.post.gz"))
		}else{
			div = fread(paste0("/home/jmurga/mkt/202004/rawData/simulations/",model,"/ne500/",nd,"/div.tsv"))
			df1 = fread(paste0(path,"/",nd,"/","out_1.0.post.gz"))
			df2 = fread(paste0(path,"/",nd,"/","out_1.0.post.gz"))
		}

		df1$alphaWeak= alphaWeak
		df1$B = bgs
		df1$analysis = n
		names(df1) = c("Posterior alpha[w]","Posterior alpha[s]","Posterior alpha","alphaWeak","B","analysis")
		tmp1 = suppressWarnings(melt(df1))
		density[[n]] = tmp1

		df2$alphaWeak= alphaWeak
		df2$B = bgs
		df2$analysis = n
		names(df2) = c("Posterior alpha[w]","Posterior alpha[s]","Posterior alpha","alphaWeak","B","analysis")
		tmp2 = suppressWarnings(melt(df2))

		rejection[[n]] = tmp2

		aw = div$dw/div$di
		as = div$ds/div$di
		a = aw + as
		
		tmpSegment = data.table(x=c(aw,as,a),xend=c(aw,as,a),y=rep(0,3),yend=rep(Inf,3),c=c("#30504f", "#e2bd9a", "#ab2710"),g=c("True alpha[w]","True alpha[s]","True alpha"))
		tmpSegment$alphaWeak = alphaWeak
		tmpSegment$B = bgs
		tmpSegment$analysis = n
		segment[[n]] = tmpSegment

		p = ggplot(tmp1) + 
			geom_density(aes(x=value,fill=variable),alpha=0.75) + 
			scale_fill_manual(values = paletteSanMiguel,labels=c(expression(paste("Posterior ",alpha[w])), expression(paste("Posterior ",alpha[s])),expression(paste("Posterior ",alpha)))) +
			geom_segment(data=tmpSegment,aes(x=x,xend=xend,y=y,yend=yend,colour=g),size=1, linetype=5) +
			scale_color_manual("True values",values=c("True alpha[w]"=paletteSanMiguel[1],"True alpha[s]"=paletteSanMiguel[2],"True alpha"=paletteSanMiguel[3]),labels=c(expression(paste("True ",alpha[w])), expression(paste("True ",alpha[s])),expression(paste("True ",alpha)))) + theme_bw() + labs(fill = "Posterior distributions",title=n,y="Posterior",x=expression(alpha))
		# ggsave(p,filename=paste0("/home/jmurga/mkt/202004/results/abc/",model,"/",n,".svg"),dpi=300)
	}

	dataPlot1 = rbindlist(density)	
	dataPlot2 = rbindlist(rejection)	
	segment = rbindlist(segment)
	
	# dataPlot1$alphaWeak = factor(dataPlot1$alphaWeak, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))
	# segment$alphaWeak = factor(segment$alphaWeak, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))
	# dataPlot1$B = factor(dataPlot1$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))
	# segment$B = factor(segment$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))

	dataPlot1$alphaWeak = as.factor(dataPlot1$alphaWeak)
	segment$alphaWeak = as.factor(segment$alphaWeak)
	segment$B  = as.factor(segment$B)
	dataPlot1$B = as.factor(dataPlot1$B)

	dataPlot2$alphaWeak = as.factor(dataPlot2$alphaWeak)
	segment$alphaWeak = as.factor(segment$alphaWeak)
	segment$B  = as.factor(segment$B)
	dataPlot2$B = as.factor(dataPlot2$B)

	if (global == TRUE){

		dataPlot1 = dataPlot1[variable == 'Posterior alpha']
		segment = segment[g == 'True alpha']

		dataPlot1$variable = as.factor(dataPlot1$variable)

		# dataM = dataM[ ,c('Mode alpha','analysis')]

		sf = scale_fill_manual(values = paletteSanMiguel[3],labels=c(expression(paste("Posterior ",alpha[w])), expression(paste("Posterior ",alpha[s])),expression(paste("Posterior ",alpha))))

		sc = scale_color_manual("True values",values=c("True alpha"=paletteSanMiguel[3]),labels=c(expression(paste("True ",alpha))))

		# od = c("analysis","True alpha","Mode alpha","Posterior alpha","alphaWeak","B")
	}else{
		sf = scale_fill_manual(values = paletteSanMiguel,labels=c(expression(paste("Posterior ",alpha[w])), expression(paste("Posterior ",alpha[s])),expression(paste("Posterior ",alpha))))

		sc = scale_color_manual("True values",values=c("True alpha[w]"=paletteSanMiguel[1],"True alpha[s]"=paletteSanMiguel[2],"True alpha"=paletteSanMiguel[3]),labels=c(expression(paste("True ",alpha[w])), expression(paste("True ",alpha[s])),expression(paste("True ",alpha)))) 

		# od = c("analysis","True alpha[w]","Mode alpha[w]","Posterior alpha[w]","True alpha[s]","Mode alpha[s]","Posterior alpha[s]","True alpha","Mode alpha","Posterior alpha","alphaWeak","B")
	}

	p2 = ggplot(dataPlot1) + geom_density(aes(x=value,fill=variable),alpha=0.75) +  sf + geom_segment(data=segment,aes(x=x,xend=xend,y=y,yend=yend,colour=g),size=1, linetype=5) + sc + facet_grid(~analysis) + theme_bw() + facet_grid(B~alphaWeak,labeller=label_parsed) + labs(fill = "Posterior distributions")

	p3 = ggplot(dataPlot2) + geom_density(aes(x=value,fill=variable),alpha=0.75) +  sf + geom_segment(data=segment,aes(x=x,xend=xend,y=y,yend=yend,colour=g),size=1, linetype=5) + sc + facet_grid(~analysis) + theme_bw() + facet_grid(B~alphaWeak,labeller=label_parsed) + labs(fill = "Posterior distributions")

	# Tables
	dataPlot1$value = round(dataPlot1$value,3)
	dataPlot2$value = round(dataPlot2$value,3)
	segment$x      = round(segment$x,3)
	
	d1 = dataPlot2 %>% group_by(analysis,variable,alphaWeak,B) %>% summarize(q=paste0(round(mean(value),3)," [",quantile(value,c(0.1)),"-",quantile(value,0.9),"]"))
	d2 = reshape2::dcast(d1,analysis~variable)
	d3 = reshape2::dcast(segment,analysis~g,value.var='x')

	d4 = merge(d2,d3,by='analysis') %>% as.data.table

	tmp = do.call(rbind,strsplit(d2$analysis,"_"))[,3:4]
	d4$alphaWeak = tmp[,1]
	d4$B = tmp[,2]
	d4 = d4 %>% arrange(B,alphaWeak)

	diff = dataPlot1 %>% group_by(analysis,variable,alphaWeak,B) %>% summarize_all(mean)
	diff2 = reshape2::dcast(diff,analysis~variable) %>% as.data.table
	diffTrue = d3[,c(1,2)]
	diff3 = cbind(as.data.table(diff2[,1]),round(abs(diff2[,2:ncol(diff2)] - diffTrue[,2:ncol(diffTrue)]),3))
	names(diff3) = c('analysis','Diff alpha[w]','Diff alpha[s]','Diff alpha')

	dout = merge(d4,diff3,by='analysis')

	dout = dout[,c("analysis","Posterior alpha[w]","True alpha[w]","Diff alpha[w]","Posterior alpha[s]","True alpha[s]","Diff alpha[s]","Posterior alpha","True alpha","Diff alpha")]


	# ggsave(p2,filename=paste0("/home/jmurga/mkt/202004/results/abc/",model,"/noDemog.svg"),dpi=300)
	out[['plot']] = p2
	out[['plotRejection']] = p3
	out[['table']] = dout
	return(out)
}

plotSimulatedAlphas = function(df,trues,title,output){

	trues = as.data.table(trues[,c('trueAlpha','path')])
	tmp  = do.call(rbind,strsplit(trues$path,"_"))
	trues$B = tmp[,4]
	trues$alphaW = tmp[,3]
	if(length(unique(df$alphaW)) > 1){
	    df$alphaW = factor(df$alphaW, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))
	    trues$alphaW = factor(trues$alphaW, labels = c("alpha[w]:0.1","alpha[w]:0.2","alpha[w]:0.3"))
	}else{
	    df$alphaW = factor(df$alphaW, labels = c("alpha[w]:0.1"))
	    trues$alphaW = factor(trues$alphaW, labels = c("alpha[w]:0.1"))	    
	}	
	if(length(unique(df$B)) > 1){
	    df$B = factor(df$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))
	    trues$B = factor(trues$B, labels = c("B:0.2","B:0.4","B:0.8","B:0.999"))    
	}else{
	    df$B = factor(df$B, labels = c("B:0.999"))
	    trues$B = factor(trues$B, labels = c("B:0.999"))
	}

	p = ggplot(df, aes(x=f,y=alphas,color=sfs)) + 
	    geom_line() + 
	    geom_hline(data=trues,aes(yintercept=trueAlpha),color = 'gray',linetype = "dotted", size = 1) +
	    facet_grid(B~alphaW,labeller=label_parsed) + 
	    ylim(-0.5,0.5) + 
	    scale_color_manual(values=c('black','#ab2710'),name = "SFS", labels = c("All alleles","Neutral + deleterious")) + 
	    theme_bw() + 
	    ylab(expression(alpha)) + 
	    xlab('Derived Allele Count') + 
	    theme(legend.position="bottom",legend.text=element_text(size=14),plot.title=element_text(hjust=0.5,face='bold'))

	ggsave(p,filename=output,dpi=600)
	return(p)
}

plotMethodsABC = function(){

	df1 = fread('rejection_1.tsv.gz');df1$analysis = 'rejection'
	df2 = fread('neuralnet_1.tsv.gz');df2$analysis = 'neuralnet'
	df3 = fread('out.0.log.post.gz');df3$analysis = 'log'
	df4 = fread('out.0.tangent.post.gz');df4$analysis = 'tangent'

	df = rbind(df1,df2,df3,df4)
	d = melt(df,id.vars='analysis')
	d$analysis = factor(d$analysis,levels=unique(d$analysis))
	p = ggplot(d) + geom_density(aes(x=value,fill=variable),alpha=0.5) + fillSanMiguel() + facet_wrap(~analysis,ncol=2)
	return(p)
}

