out_file="simu_mcexp_1"
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
include("mcexp_unmodified_functions.jl")
include("mcexp_data_initializer.jl")
include("mcexp_simulation.jl")
include("mcexp_mcmc_utils.jl")
include("mcexp_unmodified_functions.jl")
include("mcexp_wrapper.jl")
include("mcexp_burn.jl")
include("mcexp_mcmc.jl")

σ²,α, m = map(rexp, ones(3)*10)
tree_file  = out_file*".tre"
write_nexus(30, 1., 1., 0., tree_file) 

tip_values, tip_areas, tree, bts = MC(1., tree_file,σ² =  σ²,m =  0., α = 0.)

R, δt= mcexp(tip_values, tip_areas, tree, bts, out_file, 0., 0.)
