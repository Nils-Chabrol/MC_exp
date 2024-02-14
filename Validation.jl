
using Random
using Statistics
using Distributions
using Dates
using BenchmarkTools
using DelimitedFiles
using CSV
using DataFrames
using Distributed
using Plots
using LaTeXStrings
using Optim
using Tapestree
using RCall
using NewickTree
using Optim: minimizer, optimize, Options
using Random: randexp, shuffle!
using DelimitedFiles: readdlm, writedlm
using ProgressMeter: Progress, next!
using Optim: minimizer, optimize, Options
using LinearAlgebra: BLAS.axpy!, BLAS.gemm!, eigvecs, eigvals, diagm
using RCall: reval, @rput
cd("Documents/julia/MC_exp")
include("mcexp_unmodified_functions.jl")
include("mcexp_data_initializer.jl")
include("mcexp_simulation.jl")
include("mcexp_mcmc_utils.jl")
include("mcexp_unmodified_functions.jl")
include("mcexp_wrapper.jl")
include("mcexp_burn.jl")
include("mcexp_mcmc.jl")

function write_nexus(n 	      ::Int64  =  10,
					scale     ::Float64 = 1. ,
					b         ::Float64 = 1.,
					d         ::Float64 = 0.,
					file_name ::String = "test.tre")
	reval("""
    library(\"phytools\")
    tree     <- pbtree( n = $n, b = $b, d = $d, scale = $scale)
    write.nexus(tree, file= '$file_name', translate = F)
    """)
    return nothing
end

function read_nexus(tree_file      ::String; 
                    order          ::String = "cladewise", 
                    branching_times::Bool = true)

    str = reval("""
        library(\"phytools\")
        tree     <- read.nexus('$tree_file') 
        tree     <- reorder(tree, order = '$order')
        edge     <- .subset2(tree,'edge')
        Nnode    <- .subset2(tree,'Nnode')
        tiplabel <- .subset2(tree,'tip.label')
        edlength <- .subset2(tree,'edge.length')
        list(edge,Nnode,tiplabel,edlength)
        """)

    edge     = rcopy(str[1])
    edge     = convert(Array{Int64},edge)
    Nnode    = rcopy(str[2])
    Nnode    = convert(Int64,Nnode)
    tiplabel = rcopy(str[3])
    edlength = rcopy(str[4])
    edlength = convert(Array{Float64},edlength)

    tree = rtree(edge, edlength, tiplabel, Nnode)

    if branching_times
        brtimes = reval("""
            brtimes <- branching.times(tree)
            """)
        brtimes = rcopy(brtimes)
        return tree, brtimes
    else
        return tree
    end
end





function Validation(n ::Int64,
                    scale ::Float64,
                    b::Float64,
                    d::Float64,
                    tree_file::String)
	σ²,α, m = map(rexp, ones(3)*10)
	write_nexus(n, scale, b, d, tree_file) 
	tip_values, tip_areas, tree, bts = MC(1., tree_file,σ² =  σ²,m =  m, α = α)
    out_file = replace(tree_file, ".tre" =>"")
	R= mcexp(tip_values, tip_areas, tree, bts, out_file)
    return R
end


σ²,α, m = map(rexp, ones(3)*10)
tree_file  = "treest.tre"
write_nexus(30, 1., 1., 0., "treest.tre") 
tip_values, tip_areas, tree, bts = MC(1., tree_file,σ² =  σ²,m =  m, α = α)

tip_values, tip_areas, tree, bts, plotMC = plot_MC(1., tree_file, σ² =  σ²,m =  m, α = .5)
plot(plotMC)

out_file = "BM+m"
R, δt= mcexp(tip_values, tip_areas, tree, bts, out_file)

times = cumsum(δt)
pushfirst!(times, 0.)

Xc=R[6]
plotMCMC=plot(times,Xc)


plot(plotMCMC,plotMC, legend = false)

plot(plotMC)