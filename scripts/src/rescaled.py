tmpSim = pd.read_csv(PATH + "/rawData/simulations/tennesen.tsv", sep='\t')
tmpSim['path'] = tmpSim.apply(lambda row: "/home/jmurga/mkt/202004/rawData/simulations/tennesen/tennesen_" + str(row['alpha']) + "_" + str(row['alphaW']) + "_" + str(np.round(row['B'],4)),axis=1)
tmpSim['analysis'] = 'tennesen'
b = tmpSim[tmpSim.alphaW==0.3].head(2)
tmpSim = pd.concat([b,tmpSim.tail(6)])


rescaled = tmpSim.copy()
rescaled['path'] = rescaled.apply(lambda row: "/home/jmurga/mkt/202004/rawData/simulations/rescaled/tennesen_" + str(row['alpha']) + "_" + str(row['alphaW']) + "_" + str(np.round(row['B'],4)),axis=1)
rescaled['analysis'] = 'rescaled'
# rescaled = rescaled[simTable.alphaW==0.3]
simTable = pd.concat([tmpSim, rescaled])

df,dfAlpha = saveSimulatedAlphas(table=simTable);dfAlpha['B'] = round(dfAlpha['B'],3)
dfAlpha['group'] = dfAlpha['analysis'].str.capitalize() + ": " + dfAlpha.sfs


dfAlpha = dfAlpha.sort_values(['B', 'alphaW']).reset_index(drop=True)

dfAlpha.path = pd.Categorical(dfAlpha.path, categories=dfAlpha.path.unique())
p = (ggplot(dfAlpha, aes(x='f',y='alphas',color='group')) + geom_line() + facet_wrap('~path',ncol=3,scales='free_x') + scale_color_manual(['black','red','gray','#ab2710']) + theme_bw() + theme(legend_position='bottom'));p
# ggsave(p,filename="/home/jmurga/rescaledSimulations.svg")


