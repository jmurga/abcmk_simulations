using MKtest, CSV, DataFrames 

function analytical_approach(param)

	B = param.B

	x,y = MKtest.analytical_alpha(param)

	mu = param.θ_flanking/(4*param.N)

	return hcat(mu,param.ppos_l,param.ppos_h,param.al_low,param.al_tot,round(y[end],digits=5),param.B)
end

function sim_table(alphas,bgs_values,p_size,n_size,l)

	out = zeros(size(alphas,1)*size(bgs_values,1)*3,7)
	it = 1

	adap = MKtest.parameters(N=p_size,n=n_size,gam_dfe=-457,Lf=l,B=0.999,gL=10,gH=500,al_tot=0.4,al_low=0.2,B_bins=bgs_values)

	for a in 1:size(alphas,1)
		adap.al_tot = alphas[a]
		for b in 1:size(bgs_values,1)

			adap.B = bgs_values[b]

			adap.al_low = alphas[a]*0.25
			r1 = analytical_approach(adap)
			adap.al_low = alphas[a]*0.5
			r2 = analytical_approach(adap)
			adap.al_low = alphas[a]*0.75
			r3 = analytical_approach(adap)
				
			out[it,:] = r1
			out[it+1,:] = r2
			out[it+2,:] = r3
			it = it + 3
	  end
	end
	out = sort(DataFrame(out,:auto),5)
	out[:,4] = round.(out[:,4],digits=2)

	rename!(out, [:bgs_θ,:ppos_l,:ppos_h ,:alpha_weak, :alpha,:estimation,:B])
	
	return out
end

alpha = [0.4]
bgs   = [0.2,0.4,0.8,0.999]
# bgs   = [0.999] |> permutedims

simulations = sim_table(alpha,bgs,parse(Int,ARGS[1]),parse(Int,ARGS[2]),parse(Int,ARGS[3]))

CSV.write("/home/jmurga/mkt/202004/raw_data/simulations/" * ARGS[4] * ".tsv", simulations, delim='\t')

