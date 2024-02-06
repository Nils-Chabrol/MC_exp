
reps_per_period(br_length::Float64, const_δt::Float64) = 
round(Int64,cld(br_length,const_δt))


"""
    nconstant_sim!(Xt   ::Array{Float64,1},
                        nreps::Int64,
                        δt   ::Float64,
                        m    ::Float64, 
                        rate ::Float64,
                        α    ::Float64,
                        ψ    ::Float64,
                        θ    ::Float64)

Simulate the trait evolution along a speciation waiting time.
"""
function nconstant_sim!(Xt   ::Array{Float64,1},
    nreps::Int64,
    δt   ::Float64,
    m    ::Float64, 
    rate ::Float64,
    α    ::Float64,
    ψ    ::Float64,
    θ    ::Float64)

  # n species and narea areas
n = size(Xt)[1]
nch      = zeros(Int64, n)

  # allocate memory for lineage averages and differences
δX = zeros(n, n)      # X pairwise differences
la = zeros(n)         # lineage averages
for i in Base.OneTo(nreps)

    # estimate area and lineage averages and area occupancy
    δX_la!(δX, la, Xt, α, n)

    # trait step
    traitsam_1step!(Xt, la, δt, rate, m, α, ψ, θ, n)

end

return nothing
end

function nconstant_sim2!(Xt   ::Array{Float64,1},
    X_values::Vector{Vector{Float64}},  
    Time_values::Vector{Vector{Float64}},
    Time::Float64,
    nreps::Int64,
    δt   ::Float64,
    m    ::Float64, 
    rate ::Float64,
    α    ::Float64,
    ψ    ::Float64,
    θ    ::Float64)

  # n species and narea areas
  n = size(Xt)[1]
  nch      = zeros(Int64, n)

  # allocate memory for lineage averages and differences

  for i in Base.OneTo(nreps)
    δX = zeros(n, n)      # X pairwise differences
    la = zeros(n)         # lineage averages

    # estimate area and lineage averages and area occupancy
    δX_la!(δX, la, Xt, α, n)
    Time+=δt
    # trait step
    traitsam_1step2!(Xt, X_values, Time_values, Time, la, δt, rate, m, α, ψ, θ, n)
  end
  return nothing
end
"""
    δX_la!(    δX   ::Array{Float64,2},
                    la   ::Array{Float64,1},
                    Xt   ::Array{Float64,1}, 
                    α    ::Float64,
                    n    ::Int64)

Estimate lineage specific averages given sympatry configuration.
"""

function δX_la!(    δX   ::Array{Float64,2},
    la   ::Array{Float64,1},
    Xt   ::Array{Float64,1}, 
    α    ::Float64,
    n    ::Int64)
  for j = Base.OneTo(n), i = Base.OneTo(n)
    i == j && continue
      # δX differences
    δX[i,j] = Xt[i] - Xt[j]
  end
  for j = Base.OneTo(n), i = Base.OneTo(n)
    i == j && continue

    la[i] += sign(δX[i,j])*exp(-α*δX[i,j]^2)
  end
  return nothing
end




"""
    traitsam_1step!(Xt  ::Array{Float64,1}, 
                         la  ::Array{Float64,1}, 
                         δt  ::Float64, 
                         rate::Float64, 
                         m  ::Float64,
                         α ::Float64, 
                         ψ ::Float64,
                         θ ::Float64,
                         n   ::Int64)

Sample one step for trait evolution history: `X(t + δt)`.
"""

function traitsam_1step!(Xt  ::Array{Float64,1}, 
 la  ::Array{Float64,1}, 
 δt  ::Float64, 
 rate::Float64, 
 m  ::Float64,
 α ::Float64, 
 ψ ::Float64,
 θ ::Float64,
 n   ::Int64)

  @inbounds begin
    for i in Base.OneTo(n)
        Xt[i] += ψ*(θ - Xt[i]) + Eδx(la[i], m, δt) + randn()*rate

    end
  end
end

function traitsam_1step2!(Xt  ::Array{Float64,1}, 
                      X_values::Vector{Vector{Float64}},  
                      Time_values::Vector{Vector{Float64}},
                      Time::Float64,
                      la  ::Array{Float64,1}, 
                      δt  ::Float64, 
                      rate::Float64, 
                      m  ::Float64,
                      α ::Float64, 
                      ψ ::Float64,
                      θ ::Float64,
                      n   ::Int64)
  @inbounds begin
    for i in Base.OneTo(n)
        Xt[i] += ψ*(θ - Xt[i]) + Eδx(la[i], m, δt) + randn()*rate
        push!(X_values[i], Xt[i])
        push!(Time_values[i], Time)
    end
  end
end
"""
    Eδx(μ::Float64, ωx::Float64, δt::Float64)

Return the expected value given the weighted average with sympatric species.
"""
Eδx(μ::Float64, m::Float64, δt::Float64) = m * μ * δt




"""
    MC(X_initial::Float64,
                        tree_file::String;
                        m        = 1.0,
                        σ²       = 0.5,
                        α        = 0.5,
                        const_δt = 1e-4,
                        ψ        = 0.0,
                        θ        = 0.0)

Simulate the Matching-Competition model while taking into account pairwise distances between trait value, and a stabilisation process.

...
# Arguments
- `X_initial::Float64`: trait starting value.
- `tree_file::String`: full path to tree file.
- `m::Float64 = 0.0`: simulated value of ``m``, the intensity of competition relative to the stabilizing and the brownian processes.
- `σ²::Float64 = 0.5`: simulated value of ``σ^2``, the intensity of the brownian process.
- `α::Float64 = 0.0`: simulated value of ``α``, the strength of the exclusive competition between the species.
- `ω0::Float64 = 0.0`: simulated value of ``ω_0``.
- `const_δt = 1e-4`: # delta t used to approximate the simulation (lower values
  are more accurate but at a computation cost).
- `ψ::Float64 = 0.0`: simulated value of ``ψ``, the strength of the stabilizing process.
- `θ::Float64 = 0.2`: simulated value of ``θ``, the optimal value of the stabilizing process.
...
"""
function MC(X_initial::Float64,
    tree_file::String;
    σ²       = 0.5,
    m        = 1.0,
    α        = 0.5,
    const_δt = 1e-4,
    ψ        = 0.0,
    θ        = 0.0)

  # bounds checks for parameters
  0.0 >= σ² && error("σ² has to be > 0.0")
  0.0 >= m && error("m has to be > 0.0")

  tree, bts = read_nexus(tree_file)

  br = branching_times(tree)

  # sort according to branching times
  brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)

  # sort branching times
  sort!(bts, rev = true)

  # add present
  push!(bts, 0.0)

  # number of speciation events
  nbt = lastindex(bts) - 1

  # calculate speciation waiting times
  swt = Array{Float64,1}(undef,nbt)
  for i in Base.OneTo(nbt)
    swt[i] = bts[i] - bts[i+1]
  end

  # initial values
  Xt = fill(X_initial, 2)

  # start of alive
  alive = sortslices(br, dims = 1, by = x -> x[1])[1:2,2]

  nalive = lastindex(alive)

  rate = sqrt(const_δt*σ²)

  # loop through waiting times
  for j in Base.OneTo(nbt)
    nreps = reps_per_period(swt[j], const_δt)

  # simulate during the speciation waiting time
    nconstant_sim!(Xt, nreps, const_δt, m,rate, α, ψ, θ)

    if j == nbt
        break
    end

    # which lineage speciates
    wsp = brs[j,2]

    # where to insert
    wti = findfirst(isequal(wsp), alive)

    # index to insert
    idx = sort(push!(collect(1:nalive), wti))

    # insert in Xt & Yt
    Xt = Xt[idx]

    # update alive
    chs = brs[wsp .== brs[:,1],2]

    alive[wti] = chs[1]
    insert!(alive, wti+1, chs[2])

    nalive = length(alive)
  end

  tip_traits = 
  Dict(convert(Int64, alive[i]) => Xt[i]   for i = Base.OneTo(nalive))
  tip_areas = 
  Dict(convert(Int64, alive[i]) => [1]   for i = Base.OneTo(nalive))
  pop!(bts)

return tip_traits, tip_areas, tree, bts

end



function plot_MC(X_initial::Float64,
    tree_file::String;
    m        = 1.0,
    σ²       = 0.5,
    α        = 5.0,
    const_δt = 1e-4,
    ψ        = 0.0,
    θ        = 0.0)

  # bounds checks for parameters
  0.0 >= σ² && error("σ² has to be > 0.0")
  0.0 >= m && error("m has to be > 0.0")

  tree, bts = read_nexus(tree_file)

  br = branching_times(tree)
  # sort according to branching times
  brs = sortslices(br, dims = 1, by = x -> x[5], rev = true)

  # sort branching times
  sort!(bts, rev = true)

  # add present
  push!(bts, 0.0)

  # number of speciation events
  nbt = lastindex(bts) - 1

  # calculate speciation waiting times
  swt = Array{Float64,1}(undef,nbt)
  for i in Base.OneTo(nbt)
    swt[i] = bts[i] - bts[i+1]
  end
  # initial values
  Xt = fill(X_initial, 2)
  X_values= [[X_initial],[X_initial]]
  Time_values = [[0.0],[0.0]]
  Time=0.0
  # start of alive
  alive = sortslices(br, dims = 1, by = x -> x[1])[1:2,2]

  nalive = lastindex(alive)

  rate = sqrt(const_δt*σ²)
  # loop through waiting times
  for j in Base.OneTo(nbt)
    nreps = reps_per_period(swt[j], const_δt)
  # simulate during the speciation waiting time
    nconstant_sim2!(Xt, X_values, Time_values,Time, nreps, const_δt, m,rate, α, ψ, θ)
    if j == nbt
        break
    end
  # which lineage speciates
    wsp = brs[j,2]
  # where to insert
    wti = findfirst(isequal(wsp), alive) #find the position of the speciating lineage

  # index to insert
    idx = sort(push!(collect(1:nalive), wti)) #renvoie un vecteur contenant la liste des positions des espèces, en dédoublant celle de la lignée se spéciant.
  # insert in Xt & Yt
    Xt = Xt[idx] #We duplicate the value of the speciating lineage at the right position. Lineages that speciate have their two daughter lineages side by side at the same position.

  # update alive
    chs = brs[wsp .== brs[:,1],2] #chs corresponds to the numbers of the daughter nodes
    alive[wti] = chs[1] #one of the daughter nodes replace the parent node.
    insert!(alive, wti+1, chs[2]) #we insert the second daughter node at his right place.
    insert!(X_values, wti+1, [])
    insert!(Time_values, wti+1, [])
    nalive = length(alive)
    Time+=swt[j]
  end

  tip_traits = Dict(convert(Int64, alive[i]) => Xt[i]   for i = Base.OneTo(nalive))

  pop!(bts)

  return plot(Time_values,X_values)
end


