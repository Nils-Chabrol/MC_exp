

#=

Tree utilities for ESSE

Ignacio Quintero Mächler

t(-_-t)

September 19 2017

=#




"""
Immutable type of an R tree `phylo` object type.
"""
struct rtree
  ed  ::Array{Int64,2}
  el  ::Array{Float64,1}
  tlab::Array{String,1}
  nnod::Int64
end





"""
    read_tree(tree_file::String; 
              order::String = "cladewise", 
              branching_times::Bool = true)

Function to read a tree using `RCall`
to call **ape** tree reading capabilities. 
"""
function read_tree(tree_file      ::String; 
                   order          ::String = "cladewise", 
                   branching_times::Bool = true)

  str = reval("""
                library(\"ape\")
                tree     <- read.tree('$tree_file') 
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






"""
    make_ape_tree(n::Int64, 
                  λ::Float64, 
                  μ::Float64; 
                  order::String = "cladewise", 
                  branching_times::Bool = true)

Make a phylogenetic tree using `phytools` in R. 
"""
function make_ape_tree(n              ::Int64, 
                       λ              ::Float64, 
                       μ              ::Float64; 
                       order          ::String = "cladewise", 
                       branching_times::Bool   = true)

  str = reval("""
                library(ape)
                library(phytools)
                tree     <- pbtree(n = $n, b = $λ, d = $μ, extant.only = TRUE)
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





"""
    maketriads(ed::Array{Int64,2})

Make edge triads given the tree. The first number is the parent, 
the second and third the children.
"""
function maketriads(ed::Array{Int64,2}; rev::Bool = false)

  # internal nodes
  ins = unique(ed[:,1])[1:(end-1)]::Array{Int64,1}

  rev && sort!(ins, rev = true)

  ed1 = ed[:,1]
  ed2 = ed[:,2]

  trios = Array{Int64,1}[]

  # for all internal nodes
  for i in ins
    daus = findall(isequal(i), ed1)
    pushfirst!(daus, findfirst(isequal(i), ed2))
    push!(trios, daus)
  end

  return trios::Array{Array{Int64,1},1}
end





"""
    abs_time_branches(el  ::Array{Float64,1}, 
                      ed  ::Array{Int64,2},
                      ntip::Int64)

Make array with absolute initial and end time for each 
branch. Time goes backwards with the present being `0.0`.
"""
function abs_time_branches(el  ::Array{Float64,1}, 
                           ed  ::Array{Int64,2},
                           ntip::Int64)

  @inbounds begin
    # make real time edges
    elrt = zeros(size(ed))

    # edges 1
    ed1 = ed[:,1]

    for i in axes(elrt,1)
      # is not tip
      if ed[i,2] > ntip
        elrt[i,2] = elrt[findfirst(isequal(ed[i,2]), ed1)::Int64,1]
        elrt[i,1] = elrt[i,2] + el[i]
      else
        elrt[i,1] = el[i]
      end
    end

  end

  return elrt
end





"""
    brts(el  ::Array{Float64,1}, 
                ed  ::Array{Int64,2},
                ntip::Int64)

Get branching times for a tree in "cladewise" order. 
Time goes backwards with the present being `0.0`.
"""
function brts(el  ::Array{Float64,1}, 
              ed  ::Array{Int64,2},
              ntip::Int64)

  @inbounds begin
    # make real time edges

    # edges 1
    ed1 = ed[:,1]
    ed2 = ed[:,2]

    ne   = lastindex(ed1)
    nn   = zeros(ntip-1)
    intn = findall(map(x -> x > ntip, ed2))

    for i in intn
      nn[ed2[i]-ntip] = nn[ed1[i] - ntip] + el[i]
    end

    # tree height
    trh = nn[ed1[ne] - ntip] + el[ne]
    for i in Base.OneTo(ntip-1)
      nn[i] = trh - nn[i]
    end

  end

  return nn
end





"""
    tree_height(el  ::Array{Float64,1}, 
                ed  ::Array{Int64,2},
                ntip::Int64)

Estimate tree height.
"""
function tree_height(el  ::Array{Float64,1}, 
                     ed  ::Array{Int64,2},
                     ntip::Int64)

  @inbounds begin
    # tree height
    th  = 0.0
    ed2 = ed[:,2]

    da::Int64 = findfirst(isequal(1), ed2)

    # if the first branch reaches the tree height
    ed[da,1] == (ntip+1) && return el[da]

    while true
      th += el[da]

      pr = findfirst(isequal(ed[da,1]), ed2)::Int64

      if ed[pr,1] == (ntip+1)
        th += el[pr]
        break
      end

      da = findfirst(isequal(ed[pr,2]), ed2)::Int64
    end
  end

  return th
end





"""
    tip_dictionary(tS::Array{Int64,1})

Create a dictionary. WARNING: ONLY FOR USE WITHOUT CARING ABOUT TOPOLOGY.
"""
function tip_dictionary(tS::Array{Int64,1})

  # make tip values Dictionary
  tv = Dict{Int64, Int64}()

  for i in Base.OneTo(lastindex(tS))
    push!(tv, i => tS[i])
  end

  return tv
end




"""
    numberedges(ed::Array{Int64,2}, tN::Array{Int64,1})

Change numbering scheme so that tips are `1:numerofspecies` followed
by node numbering. MRCA is `numerofspecies+1`.
"""
function numberedges(ed::Array{Int64,2}, tN::Array{Int64,1}, tS::Array{Int64,1})
  nt = 1

  ni = lastindex(tN) + 1

  edc = zeros(Int64,size(ed))
  edc[1,1] = edc[2,1] = ni
  ni += 1

  e1 = ed[:,1]

  # make tip values Dictionary
  tv = Dict{Int64, Int64}()

  for i in axes(ed,1)

    nn = ed[i,2]
    ww = findfirst(isequal(nn), tN)

    # if tip
    if !isnothing(ww)
      edc[i,2] = nt
      push!(tv, nt => tS[ww])
      nt += 1
    else
      ii = findfirst(isequal(nn), e1)
      edc[ii,1] = edc[ii+1,1] = edc[i,2] = ni
      ni  += 1 
    end

  end

  return edc, tv
end




"""
    postorderedges(ed  ::Array{Int64,2},
                   el  ::Array{Float64,1},
                   ntip::Int64)

Organize edges, edge lengths in postorder traversal fashion.
"""
function postorderedges(ed  ::Array{Int64,2},
                        el  ::Array{Float64,1},
                        ntip::Int64)

  # post-order transversal using 2 stacks
  s1 = [ed[1]]
  s2 = Int64[]

  while lastindex(s1) > 0
    nod = pop!(s1)
    push!(s2, nod)

    wn = findfirst(isequal(nod), ed[:,1])
    if isnothing(wn)
      continue
    else
      push!(s1, ed[wn,2],ed[(wn+1),2])
    end
  end

  # rearrange edges accordingly
  indx = deleteat!(indexin(reverse(s2), ed[:,2]), size(ed,1)+1)

  ed = ed[indx,:]
  el = el[indx]

  # advance nodes with only daughter tips
  tnd = Int64[]
  ndp = Int64[]
  for nd in unique(ed[:,1])
    fed = findall(ed[:,1 ] .== nd)
    if length(filter(x -> x <= ntip, ed[fed,2])) == 2
      push!(tnd, nd)
      push!(ndp, fed[1], fed[2])
    end
  end

  append!(ndp,setdiff(1:(2ntip-2), ndp))

  ed[ndp,:]

  return ed[ndp,:], el[ndp]
end





"""
    remove_extinct(ed::Array{Int64,2}, 
                   el::Array{Float64,1}, 
                   ee::Array{Int64,1})

Remove extinct nodes from tree given extinct edges `ee`.
"""
function remove_extinct(ed::Array{Int64,2}, 
                        el::Array{Float64,1}, 
                        ee::Array{Int64,1})

  # identify extinct nodes
  nse = ed[ee,2]

  for i in nse

    @views e1 = ed[:,1]
    @views e2 = ed[:,2]
    r = findfirst(x -> x === i, e2)

    # parent node
    pn = e1[r] 

    # other row
    if lastindex(e1) === r
      or = r-1
    else
      or = e1[r+1] === pn ? r+1 : r-1
    end

    # ancestral row
    ar = findfirst(x -> x === pn, e2) 

    if !isnothing(ar)
      # assign node
      e2[ar]  = e2[or]
      # add edge length
      el[ar] += el[or]
    end

    # remove
    news = setdiff(1:length(e2), r, or)
    ed   = ed[news,:]   # remove from edges
    el   = el[news]     # remove from edge lengths
  end

  return ed, el
end



#=Utility functions for Tapestree

Ignacio Quintero Mächler

February 6 2017

t(-_-t)

=#


"""
    rowind(x::Int64, nrow::Int64)

Get row indexing from matrix indexing.
"""
rowind(x::Int64, nrow::Int64) = mod1(x,nrow)





"""
    colind(x::Int64, nrow::Int64)

Get column indexing from matrix indexing
"""
colind(x::Int64, nrow::Int64) = cld(x, nrow)





"""
    vecind(row::Int64, col::Int64, nrow::Int64)

Get vector indexing from column and row.
"""
vecind(row::Int64, col::Int64, nrow::Int64) = row + nrow*(col - 1)





"""
    uniroot(f, approx = 1e-8, a = 0.0, b = 0.1)

Find the root of function between `0.0` and `b`.
"""
function uniroot(f; approx = 1e-8, a = 0.0, b = 0.1) 
  # choose b
  while sign(f(a)::Float64)::Float64 == sign(f(b)::Float64)::Float64
    b += 0.1
  end
  m::Float64 = (a + b)/2.0::Float64

  while abs(f(m)::Float64)::Float64 > approx
    if sign(f(a)::Float64)::Float64 == sign(f(m)::Float64)::Float64
      a = m::Float64
    else 
      b = m::Float64
    end
    m = (a + b)/2.0::Float64
  end
  return m::Float64
end 

#=

utilities for tribe

Ignacio Quintero Mächler

t(-_-t)

May 01 2017

=#



"""
    σ²ϕprop()

Generate proportional proposals for σ² 
using random samples from **LogNormal** distributions. 
"""
σ²ϕprop() = exp(randn() - 0.1)







"""
    rexp(λ::Float64)

Generate one random sample from a **Exponential** distribution
with mean `λ`. 
"""
rexp(λ::Float64) = (randexp()/λ)::Float64






"""
    idxlessthan(x::Array{Float64,1}, val::Float64)

Get index in sorted vector `x` corresponding to the value 
that is closest to but less than `val` in sorted arrays 
using a sort of uniroot algorithm.
"""
function idxlessthan(x::Array{Float64,1}, val::Float64) 
  
  @inbounds begin

    a::Int64 = 1
    b::Int64 = lastindex(x)
  
    if x[b] < val
      return b
    end

    mid::Int64 = div(b,2)  

    while b-a > 1
      val < x[mid] ? b = mid : a = mid
      mid = div(b + a, 2)
    end

  end

  return a
end 





"""
    indmindif_sorted(x::Array{Float64,1}, val::Float64)

Get index in sorted vector `x` corresponding to the value 
that is closest to `val` in sorted arrays 
using a sort of uniroot algorithm.
"""
function indmindif_sorted(x::Array{Float64,1}, val::Float64) 
  a::Int64   = 1
  b::Int64   = lastindex(x)
  mid::Int64 = div(b,2)  

  while b-a > 1
    val < x[mid] ? b = mid : a = mid
    mid = div(b + a, 2)
  end

  abs(x[a] - val) < abs(x[b] - val) ? a : b
end 





"""
    make_edgeind(childs::Array{Int64,1}, B::Array{Float64,2})
  
Make ragged array with indexes for each edge in `Y`.
"""
function make_edgeind(childs::Array{Int64,1}, B::Array{Float64,2}, ntip::Int64)

  Bli = LinearIndices(B)
  bridx = UnitRange{Int64}[]
  for b in childs
    bindices = Bli[findall(isequal(b), B)]
    if b != (ntip+1)
      pushfirst!(bindices, bindices[1] - 1)
    end
    bidx = bindices[1]:bindices[end]
    push!(bridx, bidx)
  end

  bridx
end





"""
    make_edgeδt(bridx::Array{Array{Int64,1},1}, δt::Array{Float64,1}, m ::Int64)

Make ragged array of the cumulative δtimes for each branch.
"""
function make_edgeδt(bridx::Array{UnitRange{Int64},1}, 
                     δt   ::Array{Float64,1}, 
                     m    ::Int64)
  
  brδt = Array{Float64,1}[]
  
  for j in 1:(length(bridx)-1)
    bi = collect(bridx[j][1:(end-1)])
    for i in eachindex(bi)
      bi[i] = rowind(bi[i], m)
    end
    push!(brδt, pushfirst!(cumsum(δt[bi]),0))
  end
  
  brδt
end





"""
    create_wcol(X::Array{Float64,2})

Returns indices for columns along `m` timesteps.
"""
function create_wcol(X::Array{Float64,2})

  wNaN_x = map(!isnan, X)

  # make ragged array for non-NaN columns
  wcol = Array{Int64,1}[]
  for i = Base.OneTo(size(X,1))
    push!(wcol,findall(wNaN_x[i,:]))
  end

  return wcol
end







"""
    Eδx(μ::Float64, ωx::Float64, δt::Float64)

Return the expected value given the weighted average with sympatric species.
"""
Eδx(μ::Float64, ωx::Float64, δt::Float64) = ωx * μ * δt








"""
    makellf(δt::Vector{Float64}, Y::Array{Int64, 3}, ntip::Int64, narea::Int64)

Make likelihood and likelihood ratio functions 
for all trait matrix and biogeography history.
"""
function makellf(δt   ::Array{Float64,1}, 
                 Y    ::Array{Int64,3}, 
                 ntip ::Int64, 
                 narea::Int64,
                 m    ::Int64,
                 nedge::Int64)

  # get initial range
  wf23 = Int64[]
  for j = Base.OneTo(ntip)
    push!(wf23, findfirst(Y[:,j,1] .!= 23))
  end

  # number of normal evaluations
  n = 0
  for i = wf23
    n += length(i:(m-1))
  end

  # normal constant
  normC = -0.5*log(2.0π)*n

  function llf(X     ::Array{Float64,2},
               Y     ::Array{Int64,3}, 
               LA    ::Array{Float64,2},
               LD    ::Array{Float64,3},
               ωx    ::Float64,
               ω1    ::Float64,
               ω0    ::Float64,
               λ1    ::Float64,
               λ0    ::Float64,
               stemev::Vector{Vector{Float64}},
               brs   ::Array{Int64,3},
               σ²    ::Float64)

    ll::Float64 = normC

    @inbounds begin

      # trait likelihood
      for j = Base.OneTo(ntip)
        @simd for i = wf23[j]:(m-1)
          ll += logdnorm_tc(X[(i+1),j], 
                            X[i,j] + Eδx(LA[i,j], ωx, δt[i]), 
                            δt[i]*σ²)::Float64
        end
      end

      # biogeograhic likelihood
      for k = Base.OneTo(narea)
        ll += brll(stemev[k], λ1, λ0, brs[nedge,1,k])::Float64
        for j = Base.OneTo(ntip)
          ll += bitvectorll(Y, λ1, λ0, ω1, ω0, LD, δt, 
                            j, k, wf23[j], m)::Float64
        end
      end

    end

    return ll::Float64
  end

  function llrf(Xc     ::Array{Float64,2},
                Xp     ::Array{Float64,2},
                Yc     ::Array{Int64,3}, 
                Yp     ::Array{Int64,3}, 
                LAc    ::Array{Float64,2},
                LAp    ::Array{Float64,2},
                LDc    ::Array{Float64,3},
                LDp    ::Array{Float64,3},
                ωx     ::Float64,
                ω1     ::Float64,
                ω0     ::Float64,
                λ1     ::Float64,
                λ0     ::Float64,
                stemevc::Vector{Vector{Float64}},
                stemevp::Vector{Vector{Float64}},
                brs    ::Array{Int64,3},
                brsp   ::Array{Int64,3},
                σ²     ::Float64)

    llr::Float64 = 0.0

    @inbounds begin

      # trait likelihood
      for j = Base.OneTo(ntip)
        @simd for i = wf23[j]:(m-1)
          llr += llrdnorm_xμ(Xp[(i+1),j], Xc[(i+1),j],
                             Xp[i,j] + Eδx(LAp[i,j], ωx, δt[i]), 
                             Xc[i,j] + Eδx(LAc[i,j], ωx, δt[i]),
                             δt[i]*σ²)::Float64
        end
      end

      # biogeograhic likelihood
      for k = Base.OneTo(narea)
        llr += brll(stemevp[k], λ1, λ0, brsp[nedge,1,k]) -
               brll(stemevc[k], λ1, λ0,  brs[nedge,1,k])
        for j = Base.OneTo(ntip)
          llr += bitvectorll(Yp, λ1, λ0, ω1, ω0, LDp, δt, j, k, wf23[j], m) -
                 bitvectorll(Yc, λ1, λ0, ω1, ω0, LDc, δt, j, k, wf23[j], m)
        end
      end

    end

    return llr::Float64
  end


  return llf, llrf
end





"""
    llr_bm(Xc ::Array{Float64,2},
           Xp ::Array{Float64,2},
           δt ::Array{Float64,1},
           σ² ::Float64, 
           idx::UnitRange)

Estimate the proposal probability of a given path according to Brownian Motion.
"""
function llr_bm(Xc ::Array{Float64,2},
                Xp ::Array{Float64,2},
                idx::UnitRange,
                δt ::Array{Float64,1},
                σ² ::Float64)

  @inbounds begin
    llr::Float64 = 0.0

    @simd for i in Base.OneTo(length(idx)-1)
      llr += llrdnorm_xμ(Xc[idx[i+1]], Xp[idx[i+1]], 
                         Xc[idx[i]],   Xp[idx[i]],
                         (δt[i+1] - δt[i])*σ²)
    end
  end

  return llr::Float64
end







"""
    makellf_bgiid(bridx_a::Array{Array{UnitRange{Int64},1},1},
                  δt     ::Array{Float64,1},
                  narea  ::Int64,
                  nedge  ::Int64,
                  m      ::Int64)

Make triad likelihood function for the mutual 
independence model (iid), the proposal density 
for data augmented biogeographic histories.
"""
function makellf_bgiid(bridx_a::Array{Array{UnitRange{Int64},1},1},
                       δt     ::Array{Float64,1},
                       narea  ::Int64,
                       nedge  ::Int64,
                       m      ::Int64)

  # prepare δts
  δtA = Array{Float64,1}[]

  for j=bridx_a[1][1:(nedge-1)]
    inds = zeros(Int64,length(j) - 1)
    for i = eachindex(inds)
      inds[i] = rowind(j[i], m)
    end
    push!(δtA, δt[inds])
  end

  function fiid(Y      ::Array{Int64,3},
                stemev::Array{Array{Float64,1},1},
                brs    ::Array{Int64,3},
                triad  ::Array{Int64,1},
                λϕ1    ::Float64,
                λϕ0    ::Float64)

    ll::Float64 = 0.0

    @inbounds begin

      pr, d1, d2 = triad::Array{Int64,1}

      if pr < nedge 
        for k = Base.OneTo(narea)
          ll += bitvectorll_iid(Y, bridx_a[k][pr], λϕ1, λϕ0, δtA[pr]) +
                bitvectorll_iid(Y, bridx_a[k][d1], λϕ1, λϕ0, δtA[d1]) +
                bitvectorll_iid(Y, bridx_a[k][d2], λϕ1, λϕ0, δtA[d2])
        end
      else 
        for k = Base.OneTo(narea)
          ll += bitvectorll_iid(Y, bridx_a[k][d1], λϕ1, λϕ0, δtA[d1]) +
                bitvectorll_iid(Y, bridx_a[k][d2], λϕ1, λϕ0, δtA[d2]) +
                brll(stemev[k], λϕ1, λϕ0, brs[nedge,1,k])
        end
      end

    end

    return ll::Float64
  end


  function fiidbr(Y      ::Array{Int64,3}, 
             stemevc::Array{Array{Float64,1},1},
             brs    ::Array{Int64,3},
             br     ::Int64,
             λϕ1    ::Float64,
             λϕ0    ::Float64)

    ll::Float64 = 0.0

    @inbounds begin
      if br == nedge
        for k = Base.OneTo(narea)
          ll += brll(stemevc[k], λϕ1, λϕ0, brs[nedge,1,k])
        end
      else
        for k = Base.OneTo(narea)
          ll += bitvectorll_iid(Y, bridx_a[k][br], λϕ1, λϕ0,  δtA[br])
        end
      end
    end
    
    return ll::Float64
  end


  return fiid, fiidbr
end






"""
    stem_llr(λ1     ::Float64,
             λ0     ::Float64,
             stemc  ::Array{Int64,1},
             stemp  ::Array{Int64,1},
             stemevc::Array{Array{Float64,1},1},
             stemevp::Array{Array{Float64,1},1},
             narea  ::Int64)

Estimate likelihood ratio for stem branch.
"""
function stem_llr(λ1     ::Float64,
                  λ0     ::Float64,
                  brs    ::Array{Int64,3},
                  brsp   ::Array{Int64,3},
                  stemevc::Array{Array{Float64,1},1},
                  stemevp::Array{Array{Float64,1},1},
                  narea  ::Int64,
                  nedge  ::Int64)

  @inbounds begin

    ll::Float64 = 0.0
    for k in Base.OneTo(narea)
      ll += brll(stemevp[k], λ1, λ0, brsp[nedge,1,k]) - 
            brll(stemevc[k], λ1, λ0, brs[ nedge,1,k])
    end
  end

  return ll::Float64
end





"""
    stemiid_propr(λϕ1    ::Float64,
                  λϕ0    ::Float64,
                  stemc  ::Array{Int64,1},
                  stemp  ::Array{Int64,1},
                  stemevc::Array{Array{Float64,1},1},
                  stemevp::Array{Array{Float64,1},1},
                  narea  ::Int64) 

Proposal ratio for stem.
"""
function stemiid_propr(λϕ1    ::Float64,
                       λϕ0    ::Float64,
                       brs    ::Array{Int64,3,},
                       brsp   ::Array{Int64,3},
                       stemevc::Array{Array{Float64,1},1},
                       stemevp::Array{Array{Float64,1},1},
                       narea  ::Int64,
                       nedge  ::Int64)

  @inbounds begin

    ll::Float64 = 0.0
    for k in Base.OneTo(narea)
      ll += brll(stemevc[k], λϕ1, λϕ0, brs[ nedge,1,k]) - 
            brll(stemevp[k], λϕ1, λϕ0, brsp[nedge,1,k])
    end
  end

  return ll::Float64
end



"""
    makellr_Xupd(δt::Vector{Float64}, narea::Int64)

Make likelihood function for an internal node update in `X`.
"""
function makellr_XRupd(δt   ::Vector{Float64}, 
                       narea::Int64,
                       wcol ::Array{Array{Int64,1},1})

  δt1 = δt[1]
  wci = wcol[1]
 
  function fx(xi  ::Int64,
              xpi ::Array{Float64,1},
              X   ::Array{Float64,2},
              lapi::Array{Float64,1},
              ldpi::Array{Float64,2},
              LA  ::Array{Float64,2},
              LD  ::Array{Float64,3},
              Y   ::Array{Int64,3},
              ωx  ::Float64,
              ω1  ::Float64,
              ω0  ::Float64,
              λ1  ::Float64,
              λ0  ::Float64,
              σ²  ::Float64)

    # normal likelihoods
    llr::Float64 = 0.0

    @inbounds begin

      # loop for parent nodes
      δxim1 = δt[xi-1]
      for j = wcol[xi-1]
        llr += llrdnorm_x(xpi[j], X[xi,j], 
                          X[xi-1,j] + Eδx(LA[xi-1,j], ωx, δxim1), 
                          δxim1*σ²)
      end

      # loop for daughter nodes
      δxi = δt[xi]
      for j = wcol[xi]
        llr += llrdnorm_μ(X[xi+1, j],
                          xpi[j]  + Eδx(lapi[j],  ωx, δxi),
                          X[xi,j] + Eδx(LA[xi,j], ωx, δxi),
                          δxi*σ²)

        for k = Base.OneTo(narea)
          llr += bitbitll(Y[xi,j,k], Y[xi+1,j,k], 
                          λ1, λ0, ω1, ω0, ldpi[j,k], δxi)::Float64 -
                 bitbitll(Y[xi,j,k], Y[xi+1,j,k], 
                          λ1, λ0, ω1, ω0, LD[xi,j,k], δxi)::Float64
        end
      end
    end

    return llr::Float64
  end

  function fr(xpi ::Array{Float64,1},
             X   ::Array{Float64,2},
             σ²  ::Float64)

    llr::Float64 = 0.0

    @inbounds begin

      # loop for daughter nodes
      for j = wci
        llr += llrdnorm_μ(X[2,j], xpi[j], X[1,j], δt1*σ²)
      end
    end

    return llr::Float64
  end

  return fx, fr
end





"""
    brll(brevs::Array{Float64,1}, λ1::Float64, λ0::Float64, si::Int64)

Return likelihood for a branch in continuous time.
"""
function brll(brevs::Array{Float64,1}, λ1::Float64, λ0::Float64, si::Int64)

  ll::Float64 = 0.0

  if lastindex(brevs) > 1 
    for i = Base.OneTo(lastindex(brevs)-1)
      ll += evll(brevs[i], iszero(si) ? λ1 : λ0)::Float64
      si  = 1 - si
    end
  end

  ll += nell(brevs[end], iszero(si) ? λ1 : λ0)::Float64

  return ll::Float64
end






"""
    allλpr(λc::Array{Float64,2}, λprior::Float64)

Return log-prior for all areas 
"""
function allλpr(λ1    ::Float64,
                λ0    ::Float64,
                λprior::Float64)
  2.0*log(λprior) - λprior * (λ1 + λ0)
end

#=

data initializer for tribe model

Ignacio Quintero Mächler

t(-_-t)

June 14 2017

=#

"""
    read_data(tree_file::String, data_file::String; delim::Char = '\t', eol::Char = '\r')

Read a phylogenetic tree using **ape** package in R through 
`RCall` and the data file with the trait and biogeographic information.
"""
function read_data_tribe(tree_file::String,
                         data_file::String)

  tree, bts = read_tree(tree_file) #bts is the branching time

  tip_labels = Dict(i => val for (val,i) = enumerate(tree.tlab)) 

  data = readdlm(data_file) #renvoie un matrice any d'un jeu de donnée en fichier txt.

  if size(data,1) != (tree.nnod + 1)
    data = readdlm(data_file, '\t', '\r') 
  end

  if size(data,1) != (tree.nnod + 1)
    data = readdlm(data_file, '\t', '\n')
  end

  if size(data,1) != (tree.nnod + 1) 
    error("Data file cannot be made of the right dimensions.\n Make sure the data file has the same number of rows as tips in the tree")
  end

  data_tlab   = convert(Array{String,1}, data[:,1]) #tip labels
  data_values = convert(Array{Float64,1},data[:,2])
  data_areas  = convert(Array{Int64,2},  data[:,3:end])

  # create dictionaries
  tip_areas = Dict(tip_labels[val] => data_areas[i,:] 
                   for (i,val) = enumerate(data_tlab))

  tip_values = Dict(tip_labels[val] => data_values[i] 
                    for (i,val) = enumerate(data_tlab)) #attribute a value to a tip number

  return tip_values, tip_areas, tree, bts
end




"""
    initialize_data(tip_values::Dict{Int64,Float64}, tip_areas ::Dict{Int64,Array{Int64,1}}, m::Int64, tree::rtree, bts::Array{Float64,1})

Function to initialize `X` and `Y` matrix given
the tip_values and tip_areas (as Dictionaries).
"""
function initialize_data(tip_values::Dict{Int64,Float64},
                         tip_areas ::Dict{Int64,Array{Int64,1}},
                         min_dt    ::Float64,
                         tree      ::rtree,
                         bts       ::Array{Float64,1})

  br     = branching_times(tree) #col1 : parent node ; col2 : child node ; 3 : branch length ; 4 : parent height ; 5 : children height.
  n      = tree.nnod + 1
  nareas = length(tip_areas[1])

  #*
  # make times and δt vector
  #*

  # make sure each branch has nareas + 1 sampling points
  ets = unique(br[:,4]) #stocke les temps des époques #créer des nouvelles époques en fonction de la biogéographie.
  for i = sortperm(br[:,3])
    # number of times that cross the branch
    nover = length(findall(x -> br[i,4] > x > br[i,5], ets))
    if nover < nareas
      nets = convert(Array{Float64,1},
        range(br[i,4], stop = br[i,5], length = nareas-nover+2))
      if length(nets) > 2
        append!(ets,nets[2:(end-1)])
      end
    end
  end

  # epoch times
  sort!(ets, rev = true) #range les époques dans l'ordre.

  tr_height = ets[1]
  mδt       = min_dt*tr_height #split time into min_dt segments.

  # incorporate more 'ets' according to min_dt
  new_ets = Float64[]
  for i in eachindex(ets)

    if i == lastindex(ets)
      if ets[i]/tr_height > min_dt    
        append!(new_ets, collect(0:mδt:ets[i])[2:end])
      end
    else
      if (ets[i] - ets[i+1])/tr_height > min_dt    
        append!(new_ets, collect(ets[i+1]:mδt:ets[i])[2:end])
      end
    end
  end

  # add new_ets
  append!(ets, new_ets)

  # sort epoch times from start to end
  sort!(ets, rev = true)

  # push present
  push!(ets, 0.0)

  #create δt vector
  δt = abs.(diff(ets))

  # initialize data augmentation matrices
  X = fill(NaN, length(ets), n) #each lign correspond to the trait value at one point in time, each column represents a point in time.
  B = copy(X)
  Y = fill(23, length(ets), n, nareas)

  # which rows are branching points
  wch = indexin(bts, ets) #trouve les index des points de branchements dans les points de temps.

  # coupled nodes (cells that are coupled in array)
  coup = zeros(Int64, tree.nnod, 3)

  bord = sortperm(br[:,5])[(1:tree.nnod-1) .+ n]

  alive = collect(1:n)

  #which column alive
  wca = 1:n 

  wts = sort(wch, rev = true) .+ 1 # which time split

  setindex!(coup, wts .- 1, :,3) #place dans coup les valeurs de timesplit

  wrtf = wts[1]:length(ets) #correspond aux derniers segments de temps avant le présent (après le dernier split).
  X[wrtf, wca] .= 1.0 
  B[wrtf, wca]  = repeat(alive, inner = length(wrtf)) #remplie les même ligne que X au dessus, mais avec le numéro de la colonne (ex  pour une ligne : [1.0, 2.0, 3.0, ..., 14.0)

  for i = Base.OneTo(tree.nnod-1) #retrouve les noeuds créer à t0
    fn  = br[bord[i],2]
    wda = br[findall(isequal(fn), br[:,1]),2]

    cda = indexin(wda,alive)
    coup[i,1:2] = cda

    alive[cda] = [fn,0]
    wrtf = wts[i+1]:(wts[i]-1)

    B[wrtf,:] = repeat(alive, inner = length(wrtf))

    wca = findall(x -> x > 0, alive)

    X[wts[i+1]:length(ets),wca] .= 1.0
  end

  coup[tree.nnod,1:2] = findall(x -> x > 0, alive)
 
  ncoup = zeros(Int64,tree.nnod,2)

  for j = Base.OneTo(tree.nnod), i=1:2
    ncoup[j,i] = vecind(coup[j,3], coup[j,i], lastindex(ets))
  end

  X[1,1]      = 1.0
  B[1,1]      = n + 1
  B[B .== 0] .= NaN

  # Brownian bridges initialization for X
  si = initialize_X!(tip_values, X, B, ncoup, δt, tree)

  # declare non-23s for Y
  initialize_Y!(tip_areas, Y, B)

  return X, Y, B, ncoup, δt, tree, si
end






"""
    initialize_Y!(tip_areas::Dict{Int64,Array{Int64,1}}, Y::Array{Int64,3}, B::Array{Float64,2})

Simple function to initialize Y with all 1s (the
real biogeographic sampling is done during the MCMC).
"""
function initialize_Y!(tip_areas::Dict{Int64,Array{Int64,1}},
                       Y        ::Array{Int64,3},
                       B        ::Array{Float64,2})

  # index non 23 for Y
  ind = findall(!isnan, B)
  lB  = length(ind)

  # #D indices
  indY = CartesianIndex{3}[]
  for a = 1:size(Y,3), i = ind
    push!(indY, CartesianIndex(i,a))
  end

  Y[indY] .= 1

  for i in Base.OneTo(size(Y,2))
    Y[end,i,:] = tip_areas[i]
  end

end





"""
Maximum likelihood Brownian Motion.
"""
function make_bm_ll(tip_values::Dict{Int64,Float64},
                    tree      ::rtree)

  wt = tree.ed .<= (tree.nnod + 1) #on toruve les emplacement des espèces actuelles dans ed.
  tips_ind = findall(wt)

  # base with trait values
  ntr  = zeros(size(tree.ed))

  # assign tip values
  for i = eachindex(tip_values)
    ntr[tips_ind[i]] = tip_values[i]
  end

  # internal nodes
  ins  = unique(tree.ed[:,1])
  lins = length(ins)

  # make triads for all internal nodes
  # including the root
  trios = Array{Int64,1}[]
  ndi  = ins[1]
  daus = findall(isequal(ndi), tree.ed[:,1]) #on trouve les enfants de la racine
  pushfirst!(daus, 0)
  push!(trios, daus) #

  # for all internal nodes
  for i = 2:lins
    ndi  = ins[i]
    daus = findall(isequal(ndi), tree.ed[:,1])
    pushfirst!(daus, findall(isequal(ndi), tree.ed[:,2])[1])
    push!(trios, daus)
  end

  el = tree.el

  function f(p::Array{Float64,1})

    @inbounds begin

      σ² = p[1]

      if σ² <= 0
        return Inf
      end

      for i = eachindex(trios)
        pr, d1, d2 = trios[i] #on assigne le numéro des noeuds
        if pr == 0
          ntr[d1,1] = ntr[d2,1] = p[i+1]
        else
          ntr[pr,2] = ntr[d1,1] = ntr[d2,1] = p[i+1]
        end
      end
      
      ll = 0.0
      for i in eachindex(el)
        ll -= logdnorm_tc(ntr[i,2], ntr[i,1], el[i]*σ²)
      end
    end  
    
    return ll
  end

  return f
end





"""
    branching_times(tree::rtree)

Function to estimate absolute
branching times, with time 0 at the
present, time going backward.
"""
function branching_times(tree::rtree)

  @inbounds begin

    n    = tree.nnod + 1
    el_t = findall(x -> x <= n, tree.ed[:,2])

    brs = zeros(lastindex(tree.el),5)

    brs[:, 1:2] = tree.ed
    brs[:, 3]   = tree.el

    for j = eachindex(tree.el)
      if brs[j,1] == (n+1)
        brs[j,4] = 0.0
        brs[j,5] = brs[j,4] + brs[j,3]
      else
        brs[j,4] = brs[brs[j,1] .== brs[:,2],5][1]
        brs[j,5] = brs[j,4] + brs[j,3]
      end
    end

     # change time forward order
    @views tree_depth = brs[n .== brs[:,2],5][1]

    for j = eachindex(tree.el) 
      brs[j,4] = tree_depth - brs[j,4]
      brs[j,5] = tree_depth - brs[j,5]
    end

    brs[el_t,5] .= 0.0

  end

  return brs
end




"""
    bb(xs::Float64, xf::Float64, δt::Array{Float64,1})

Brownian bridge simulation function for
a vector of times δt.
"""
function bb(xs::Float64, xf::Float64, δt::Array{Float64,1})

  t  = pushfirst!(cumsum(δt),0.0)
  tl = lastindex(t)
  w  = zeros(tl)

  for i = Base.OneTo(tl-1)
    w[i+1] = randn()*sqrt(δt[i])
  end

  cumsum!(w, w)
  wf = w[tl]
  tf = t[tl]

  return @. xs + w - t/tf * (wf - xf + xs)
end

#=

Denisity functions for Tapestree

Ignacio Quintero Mächler

t(-_-t)

September 23 2017

=#





"""
    logdexp(x::Float64, λ::Float64)

Compute the logarithmic transformation of the 
**Exponential** density with mean `λ` for `x`.
"""
logdexp(x::Float64, λ::Float64) = (log(λ) - λ * x)::Float64




"""
    logdunifU(x::Float64)

Standard uniform distribution (`[0.0,1.0]`).
"""
logdunifU(x::Float64) = 0.0





"""
    logdunif(x::Float64, a::Float64, b::Float64)

Make function to compute the logarithmic transformation of the 
**Uniform** density with lower bound `a` and upper bound `b` for `x`.
"""
function logdunif(x::Float64, a::Float64, b::Float64)
  if x < a 
    return -Inf
  elseif x <= b 
    return -log(b-a)::Float64
  else 
    return -Inf
  end
end





"""
    llrdexp_x(xp::Float64, xc::Float64, λ::Float64)

Compute the loglik ratio of the 
**Exponential** density for proposal 
`xp` given current `xc` both with mean `λ`.
"""
llrdexp_x(xp::Float64, xc::Float64, λ::Float64) = 
  λ * (xc - xp)





"""
    logdbeta(x::Float64, α::Float64, β::Float64)

Compute the logarithmic transformation of the 
**Beta** density with shape `α` and shape `β` for `x`.
"""
logdbeta(x::Float64, α::Float64, β::Float64) = 
  ((α-1.0) * log(x)                 +
  (β-1.0) * log(1.0 - x)           +
  log(gamma(α + β)/(gamma(α)*gamma(β))))





"""
    llrdbeta_x(xp::Float64, xc::Float64, α::Float64, β::Float64)

Compute the logarithmic ratio for the **Beta** density 
with shape `α` and shape `β` between `xp` and `xc`.
"""
function llrdbeta_x(xp::Float64, xc::Float64, α::Float64, β::Float64) 
  if !(0.0 < xp < 1.0)
    return -Inf
  else
    return ((α-1.0) * log(xp/xc) +
            (β-1.0) * log((1.0 - xp)/(1.0 - xc)))
  end
end




"""
    logdnorm(x::Float64, μ::Float64, σ²::Float64)
  
Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and variance `σ²` for `x`.
"""
logdnorm(x::Float64, μ::Float64, σ²::Float64) = 
  -(0.5*log(2.0π) + 0.5*log(σ²) + (x - μ)^2/(2.0σ²))





"""
    logdnorm_tc(x::Float64, μ::Float64, σ²::Float64)

Compute the logarithmic transformation of the 
**Normal** density with mean `μ` and variance `σ²` for `x`, up to a constant
"""
logdnorm_tc(x::Float64, μ::Float64, σ²::Float64) =
  -0.5*log(σ²) - (x - μ)^2/(2.0σ²)::Float64








"""
    llrdnorm_σ²(x::Float64, μ::Float64, σ²p::Float64, σ²c::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `σ²` updates
"""
llrdnorm_σ²(x::Float64, μ::Float64, σ²p::Float64, σ²c::Float64) = #is it normal ?
  -0.5*(log(σ²p/σ²c) + (x - μ)^2*(1.0/σ²p - 1.0/σ²c))





"""
    llrdnorm_μ(x::Float64, μp::Float64, μc::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `μ` updates
"""
llrdnorm_μ(x::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  ((x - μc)^2 - (x - μp)^2)/(2.0σ²)





"""
    llrdnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `x` updates
"""
llrdnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64) =
  ((xc - μ)^2 - (xp - μ)^2)/(2.0σ²)




"""
    llrdnorm_xμ(xp::Float64, xc::Float64, μp::Float64, μc::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Normal** density 
for `x` and `μ` updates
"""
llrdnorm_xμ(xp::Float64, xc::Float64, μp::Float64, μc::Float64, σ²::Float64) =
  ((xc - μc)^2 - (xp - μp)^2)/(2.0σ²)





"""
    logdtnorm(x::Float64, σ²::Float64)

Compute the log-likelihood density for the **Truncated Normal** density
with `a = -1` and `b = Inf`
"""
function logdtnorm(x::Float64, σ²::Float64)
  if x < -1.0 
    return -Inf
  else
    return (-x^2/(2.0*σ²) - 0.5*log(2.0π) - log(0.5*sqrt(σ²)) -
            log(1.0 - erf_custom(-1.0/(sqrt(2.0σ²)))))
  end
end





"""
    llrtnorm_x(xp::Float64, xc::Float64, μ::Float64, σ²::Float64)

Compute the log-likelihood ratio for the **Truncated Normal** density
for `x` with `a = -1`, `b = Inf`
"""
function llrdtnorm_x(xp::Float64, xc::Float64, σ²::Float64)
  if xp < -1.0 
    return -Inf
  else
    return (xc^2 - xp^2)/(2.0σ²)
  end
end






"""
    erf(x::Float64)

Compute the error function
"""
function erf_custom(x::Float64)
  z = abs(x)
  t = 1.0/(1.0+0.5*z)
  r = t*exp(-z*z-1.26551223 +
        t*(1.00002368 +
            t*(0.37409196 + 
              t*(0.09678418 +
                t*(-0.18628806 +
                  t*(0.27886807 +
                    t*(-1.13520398 +
                      t*(1.48851587 +
                        t*(-0.82215223 + 
                          t*0.17087277)))))))))

  return x >= 0.0 ? (1.0-r) : (r -1.0)
end





"""
    logdhcau(x::Float64, scl::Float64)

Compute the logarithmic transformation of the 
**Half-Cauchy** density with scale `scl` for `x`.
"""
logdhcau(x::Float64, scl::Float64) = 
  log(2.0 * scl/(π *(x * x + scl * scl)))





"""
    logdhcau1(x::Float64)
  
Compute the logarithmic transformation of the 
**Half-Cauchy** density with scale of 1 for `x`.
"""
logdhcau1(x::Float64) = 
  log(2.0/(π * (x * x + 1.)))

#=

functions to estimate area and lineage trait averages 
and lineage specific differences

Ignacio Quintero Mächler

t(-_-t)

May 15 2017

=#





"""
    deltaX!(δX   ::Array{Float64,3}, 
            X    ::Array{Float64,2},
            m    ::Int64,
            ntip ::Int64,
            narea::Int64)

Estimate pairwise trait and range distances between all lineages
conditional on a given `δY`.
"""
function deltaX!(δX   ::Array{Float64,3}, 
                 X    ::Array{Float64,2},
                 wcol ::Array{Array{Int64,1},1},
                 m    ::Int64,
                 ntip ::Int64,
                 narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m), l = wcol[i], j = wcol[i]
      l == j && continue
      δX[j,l,i] = X[i,j] - X[i,l]
    end

  end

  return nothing
end





"""
    deltaY!(δY   ::Array{Float64,3}, 
            Y    ::Array{Int64,3},
            m    ::Int64,
            ntip ::Int64,
            narea::Int64)

Estimate pairwise trait and range distances between all lineages
conditional on a given `δY`.
"""
function deltaY!(δY   ::Array{Float64,3}, 
                 Y    ::Array{Int64,3},
                 wcol ::Array{Array{Int64,1},1},
                 m    ::Int64,
                 ntip ::Int64,
                 narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m), l = wcol[i], j = wcol[i]

      j == l && continue

      sl        = 0.0
      δY[l,j,i] = 0.0
      @simd for k = Base.OneTo(narea)
        if Y[i,j,k] == 1
          sl        += 1.0
          δY[l,j,i] += Float64(Y[i,l,k])
        end
      end
      δY[l,j,i] /= sl
    end

  end

  return nothing
end





"""
    deltaXY!(δX   ::Array{Float64,3}, 
             δY   ::Array{Float64,3},
             X    ::Array{Float64,2},
             Y    ::Array{Int64,3},
             m    ::Int64,
             ntip ::Int64,
             narea::Int64)

Estimate pairwise trait and range distances between all lineages.
"""
function deltaXY!(δX   ::Array{Float64,3}, 
                  δY   ::Array{Float64,3},
                  X    ::Array{Float64,2},
                  Y    ::Array{Int64,3},
                  wcol ::Array{Array{Int64,1},1},
                  m    ::Int64,
                  ntip ::Int64,
                  narea::Int64)

  @inbounds begin

    for i = Base.OneTo(m), l = wcol[i], j = wcol[i]

      j == l && continue

      # X differences
      δX[j,l,i] = X[i,j] - X[i,l]

      # Y area overlap
      sl        = 0.0
      δY[l,j,i] = 0.0
      @simd for k = Base.OneTo(narea)
        if Y[i,j,k] == 1
          sl        += 1.0
          δY[l,j,i] += Float64(Y[i,l,k])
        end
      end
      δY[l,j,i] /= sl
    end

  end

  return nothing
end






#=

Proposal functions for joint
Biogeographic competition model

Ignacio Quintero Mächler

t(-_-t)

May 16 2017

=#





#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
`Y` IID proposal functions
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#





"""
    upnode!(λ1     ::Float64,
            λ0     ::Float64,
            triad  ::Array{Int64,1},
            Y      ::Array{Int64,3},
            stemevs::Array{Array{Float64,1},1},
            bridx_a::Array{Array{UnitRange{Int64},1},1},
            brδt   ::Vector{Vector{Float64}},
            brl    ::Vector{Float64},
            brs    ::Array{Int64,3},
            narea  ::Int64,
            nedge  ::Int64)

Update node and incident branches using discrete 
Data Augmentation for all areas using a non-competitive 
mutual-independence Markov model.
"""
function upnode!(λ1     ::Float64,
                 λ0     ::Float64,
                 triad  ::Array{Int64,1},
                 Y      ::Array{Int64,3},
                 stemevs::Array{Array{Float64,1},1},
                 bridx_a::Array{Array{UnitRange{Int64},1},1},
                 brδt   ::Vector{Vector{Float64}},
                 brl    ::Vector{Float64},
                 brs    ::Array{Int64,3},
                 narea  ::Int64,
                 nedge  ::Int64)

  @inbounds begin
   
    # define branch triad
    pr, d1, d2 = triad

    # sample
    samplenode!(λ1, λ0, pr, d1, d2, brs, brl, narea)

    # save extinct
    ntries = 1
    while iszero(sum(view(brs,pr,2,:)))
      samplenode!(λ1, λ0, pr, d1, d2, brs, brl, narea)

      ntries += 1
      if ntries == 500
        return false
      end
    end

    # set new node in Y
    @simd for k in Base.OneTo(narea)
      Y[bridx_a[k][d1][1]] = Y[bridx_a[k][d2][1]] = brs[pr,2,k]
    end

    # sample a consistent history
    createhists!(λ1, λ0, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                 stemevs, brl[nedge])

    ntries = 1
    while ifextY(Y, stemevs, triad, brs, brl[nedge], narea, bridx_a, nedge)
      createhists!(λ1, λ0, Y, pr, d1, d2, brs, brδt, bridx_a, narea, nedge,
                   stemevs, brl[nedge])

      ntries += 1
      if ntries == 500
        return false
      end
    end

  end

  return true
end






"""
    samplenode!(λ1   ::Float64, 
                λ0   ::Float64,
                pr   ::Int64,
                d1   ::Int64,
                d2   ::Int64,
                brs  ::Array{Int64,3},
                brl  ::Array{Float64,1},
                narea::Int64)

Sample one internal node according to 
mutual-independence model transition probabilities.
"""
function samplenode!(λ1   ::Float64, 
                     λ0   ::Float64,
                     pr   ::Int64,
                     d1   ::Int64,
                     d2   ::Int64,
                     brs  ::Array{Int64,3},
                     brl  ::Array{Float64,1},
                     narea::Int64)
  @inbounds begin

    # estimate transition probabilities
    pr0_1, pr0_2 = Ptrfast_start(λ1, λ0, brl[pr], Val{0})
    pr1_1, pr1_2 = Ptrfast_start(λ1, λ0, brl[pr], Val{1})
    d10_1, d10_2 = Ptrfast_end(  λ1, λ0, brl[d1], Val{0})
    d11_1, d11_2 = Ptrfast_end(  λ1, λ0, brl[d1], Val{1})
    d20_1, d20_2 = Ptrfast_end(  λ1, λ0, brl[d2], Val{0})
    d21_1, d21_2 = Ptrfast_end(  λ1, λ0, brl[d2], Val{1})

    for k = Base.OneTo(narea)

      if iszero(brs[pr,1,k])
        ppr_1, ppr_2 = pr0_1, pr0_2
      else 
        ppr_1, ppr_2 = pr1_1, pr1_2
      end

      if iszero(brs[d1,2,k])
        pd1_1, pd1_2 = d10_1, d10_2
      else 
        pd1_1, pd1_2 = d11_1, d11_2
      end

      if iszero(brs[d2,2,k])
        pd2_1, pd2_2 = d20_1, d20_2
      else 
        pd2_1, pd2_2 = d21_1, d21_2
      end

      tp = normlize(*(ppr_1, pd1_1, pd2_1),
                    *(ppr_2, pd1_2, pd2_2))::Float64

      # sample the node's character
      brs[pr,2,k] = brs[d1,1,k] = brs[d2,1,k] = coinsamp(tp)::Int64
    end
  end

  return nothing
end









#=
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
X proposal functions
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
=#




"""
    uptrioX!(pr   ::Int64, 
             d1   ::Int64,
             d2   ::Int64, 
             X    ::Array{Float64,2}, 
             bridx::Array{UnitRange{Int64},1},
             brδt ::Array{Array{Float64,1},1}, 
             σ²c  ::Float64)

Update the node and adjoining branches of `trio` using Brownian bridges.
"""
function uptrioX!(pr   ::Int64, 
                  d1   ::Int64,
                  d2   ::Int64, 
                  X    ::Array{Float64,2}, 
                  bridx::Array{UnitRange{Int64},1},
                  brδt ::Array{Array{Float64,1},1},
                  brl  ::Array{Float64,1},
                  σ²ϕ  ::Float64, 
                  nedge::Int64)
  @inbounds begin

    # if not root
    if pr != nedge

      ipr = bridx[pr]
      id1 = bridx[d1]
      id2 = bridx[d2]

      # update node
      X[id1[1]] = 
      X[id2[1]] = trioupd(X[ipr[1]], 
                          X[id1[end]], 
                          X[id2[end]],
                          brl[pr], brl[d1], brl[d2], σ²ϕ)

      #update branches
      bbX!(X, ipr, brδt[pr], σ²ϕ)
      bbX!(X, id1, brδt[d1], σ²ϕ)
      bbX!(X, id2, brδt[d2], σ²ϕ)

    else
      id1 = bridx[d1]
      id2 = bridx[d2]

      # update node
      X[id1[1]] = 
      X[id2[1]] = duoupd(X[id1[end]],
                         X[id2[end]], 
                         brl[d1], brl[d2], σ²ϕ)

      # update branches
      bbX!(X, id1, brδt[d1], σ²ϕ)
      bbX!(X, id2, brδt[d2], σ²ϕ)
    end

  end

  return nothing
end





"""
    upbranchX!(j    ::Int64, 
               X    ::Array{Float64,2}, 
               bridx::Array{UnitRange{Int64},1},
               brδt ::Array{Array{Float64,1},1}, 
               σ²c  ::Float64)

Update a branch j in X using a Brownian bridge.
"""
function upbranchX!(j    ::Int64, 
                    X    ::Array{Float64,2}, 
                    bridx::Array{UnitRange{Int64},1},
                    brδt ::Array{Array{Float64,1},1},
                    σ²ϕ  ::Float64)

  @inbounds bbX!(X, bridx[j], brδt[j], σ²ϕ)

  return nothing
end





"""
    bbX!(X::Array{Float64,2}, idx::UnitRange, t::Array{Float64,1}, σ::Float64)

Brownian bridge simulation function for updating a branch in X in place.
"""
function bbX!(X  ::Array{Float64,2}, 
              idx::UnitRange,
              t  ::Array{Float64,1},
              σ²ϕ::Float64)

  @inbounds begin

    xf::Float64 = X[idx[end]]

    for i = Base.OneTo(lastindex(t)-1)
      X[idx[i+1]] = (X[idx[i]] + randn()*sqrt((t[i+1] - t[i])*σ²ϕ))::Float64
    end

    invte::Float64 = 1.0/t[end]
    xdif ::Float64 = (X[idx[end]] - xf)

    @simd for i = Base.OneTo(lastindex(t))
      X[idx[i]] = (X[idx[i]] - t[i] * invte * xdif)::Float64
    end
  end

  return nothing
end

#=

MCMC utility functions for Tapestree

Ignacio Quintero Mächler

February 6 2017

t(-_-t)

=#



"""
    uniupt(p::Float64, tn::Float64)

Uniform parameter window move.
"""
uniupt(p::Float64, tn::Float64) = abs(p + (rand()-0.5) * tn)




"""
    addupt(p::Float64, tn::Float64)

Gaussian parameter window move.
"""
addupt(p::Float64, tn::Float64) = p + randn() * tn




"""
    addupt_lims(p::Float64, tn::Float64, xmin::Float64, xmax::Float64)

Gaussian parameter window move within a region of interest using rejection.
"""
function addupt_lims(p::Float64, tn::Float64, xmin::Float64, xmax::Float64)

  s = p + randn() * tn
  if !(xmin < s < xmax)
    s = p
  end

  return s
end




"""
    addupt(p::Float64, tn::Float64)

Gaussian parameter window move to vector.
"""
function addupt!(p::Vector{Float64}, tn::Vector{Float64}, j::Int64, i::Int64)

  @inbounds p[j] += randn() * tn[i]

  return nothing
end




"""
    duoupd(xd1::Float64,
           xd2::Float64,
           td1::Float64, 
           td2::Float64,
           σ²ϕ::Float64)

Duo of Gaussians parameter update.
"""
function duoupd(xd1::Float64,
                xd2::Float64,
                td1::Float64, 
                td2::Float64,
                σ²ϕ::Float64)
  invt = 1.0/(td1 + td2)
  return randn()*sqrt(td1 * td2 * invt * σ²ϕ) + 
    (td2 * invt * xd1 + td1 * invt * xd2)
end






"""
    trioupd(xpr::Float64,
            xd1::Float64,
            xd2::Float64,
            tpr::Float64, 
            td1::Float64, 
            td2::Float64,
            σ²ϕ::Float64)

Trio of Gaussians parameter update.
"""
function trioupd(xpr::Float64,
                 xd1::Float64,
                 xd2::Float64,
                 tpr::Float64, 
                 td1::Float64, 
                 td2::Float64,
                 σ²ϕ::Float64)

    t = 1.0/(1.0/tpr + 1.0/td1 + 1.0/td2)
    return randn()*sqrt(t*σ²ϕ) + (xpr/tpr + xd1/td1 + xd2/td2)*t
end






"""
    absaddupt(p::Float64, tn::Float64)

Non-negative Gaussian parameter window move.
"""
absaddupt(p::Float64, tn::Float64) = abs(p + randn() * tn)




"""
    mulupt(p::Float64, tn::Float64)

Multiplicative parameter window move.
"""
mulupt(p::Float64, tn::Float64) = p * exp((rand() - 0.5) * tn)





"""
    makescalef(obj_ar::Float64)

Make scaling function given the objective acceptance rates.
"""
function makescalef(obj_ar::Float64)
  nar::Float64 = 1.0 - obj_ar

  function f(window::Float64, rate::Float64)
    if rate > obj_ar
      window *= (1.0 + (rate - obj_ar) / nar)::Float64
    else
      window /= (2.0 - rate / obj_ar)::Float64
    end
    
    return window::Float64
  end

  return f
end




"""
    globalscalef(λ::Float64, grate::Float64, stepsize::Float64, obj_ar::Float64)

Estimate global scaling factor.
"""
function globalscalef(λ       ::Float64, 
                      grate   ::Float64, 
                      stepsize::Float64, 
                      obj_ar  ::Float64)
  return (exp(log(λ) + stepsize * (grate - obj_ar)))::Float64
end





"""
    adaptiveupd!(Σ::Array{Float64,2}, psam::Array{Float64,1}, pmean::Array{Float64,1}, stepsize::Float64)

Adaptive update for parameter mean and Σ in place.
"""
function adaptiveupd!(Σ       ::Array{Float64,2},
                      psam    ::Array{Float64,1},
                      pmean   ::Array{Float64,1},
                      stepsize::Float64)

  @inbounds begin
    for i in Base.OneTo(length(pmean))
      psam[i]  -= pmean[i]::Float64
      pmean[i] += (stepsize * psam[i])::Float64
    end

    BLAS.axpy!(stepsize,
               BLAS.gemm!('N', 'T', 1.0, psam, psam, -1.0, copy(Σ)),
               Σ)
  end
end





"""
    makestepsize(C::Float64, η::Float64)

Make function for the stepsize for the adaptive update.
"""
function makestepsize(η::Float64)
  
  β::Float64 = rand(range((1.0/(1.0 + η)),1))

  function f(t::Float64, C::Float64)
    return (C/(t^β))::Float64
  end

  return f
end




"""
    makemvnproposal(Σ::Array{Float64,2})

Make the multivariate update given the covariance matrix.
"""
function makemvnproposal(Σ::Array{Float64,2})

  spde = *(eigvecs(Σ),sqrt.(diagm(eigvals(Σ))))
  ln   = size(spde,1)

  function f(pvec::Array{Float64,1})
    (pvec .+ 
     BLAS.gemv('N', spde, randn(ln))::Array{Float64,1}
     )::Array{Float64,1}
  end

  return f
end





















