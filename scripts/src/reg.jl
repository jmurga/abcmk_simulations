using Analytical

p, r = Analytical.ABCreg(data=ARGS[1],prior=ARGS[2], nparams=3, nsummaries=104, outputPath="/home/jmurga/dataAbc/", outputPrefix="slim1", tolerance=0.001, regressionMode="T",regPath="/home/jmurga/ABCreg/src/reg");
