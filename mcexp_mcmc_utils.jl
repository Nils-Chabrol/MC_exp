"""
    make_mhr_upd_X(Xnc1     ::Array{Int64,1},
                   Xnc2     ::Array{Int64,1},
                   wcol     ::Array{Array{Int64,1},1},
                   ptn      ::Array{Float64,1},
                   wXp      ::Array{Int64,1},
                   narea    ::Int64,
                   ntip     ::Int64,
                   Xupd_llr ::Function,
                   Rupd_llr ::Function)

Make DA update X.
"""
function make_mhr_upd_X2(X        ::Array{Float64,2},
                        Xnc1     ::Array{Int64,1},
                        Xnc2     ::Array{Int64,1},
                        wcol     ::Array{Array{Int64,1},1},
                        ptn      ::Array{Float64,1},
                        wXp      ::Array{Int64,1},
                        nstep        ::Int64,
                        δt        ::Array{Float64,1},
                        ntip     ::Int64,
                        Xupd_llr ::Function,
                        Rupd_llr ::Function)

  #Cartesian Indices
  Xcidx = CartesianIndices(X)
  rj = Xcidx[Xnc2[findfirst(isone, Xnc1)]][2] #position de la première spéciation (?)

  xpi  = fill(NaN, ntip)
  δxi  = fill(NaN, ntip, ntip)
  lani = fill(NaN, ntip)

  function f(up  ::Int64,
             Xc  ::Array{Float64,2},
             δXc ::Array{Float64,3},
             δYc ::Array{Float64,3},
             σ²c ::Float64,
             mc  ::Float64,
             αc  ::Float64,
             llc ::Float64,
             LAnc::Array{Float64,2})

    @inbounds begin

      upx = wXp[up - 3]::Int64                 # X indexing

      # if root
      if upx == 1
        # allocate
        @simd for i = Base.OneTo(ntip)
          xpi[i] =  Xc[1,i]
        end

        # update xi
        addupt!(xpi, ptn, 1, up)
        xpi[rj] = xpi[1]::Float64
        llr = Rupd_llr(xpi, Xc, σ²c)::Float64

        if -randexp() < llr #Si proposition est acceptée.
          llc    += llr::Float64
          Xc[1,:] = xpi::Array{Float64,1} #On modifie la valeur à la racine.
        end

      else

        xi, xj = Xcidx[upx].I

        # allocate
        for j = Base.OneTo(ntip)
          xpi[j]  = Xc[xi,j]
          lani[j] = LAnc[xi,j]
          @simd for i = Base.OneTo(ntip)
            δxi[i,j] = δXc[i,j,xi]
          end
        end
        xppi=Xc[xi-1, xj]
        # update xi
        addupt!(xpi, ptn, xj, up) 
        # xpi[xj] = rand(Normal(xppi+Eδx(LAnc[xi-1,xj], m, δt[xi]), δt[xi]σ²c))
          # addupt2!(xpi,xj, σ²c)
        if in(upx, Xnc1)        # if an internal node
          xpi[Xcidx[Xnc2[findfirst(isequal(upx),Xnc1)]][2]] = xpi[xj] #???
        end

        # calculate new averages
        Xupd_linavg2!(δxi, δYc, lani, wcol, xpi, xi, xj, mc, αc)
        llr = Xupd_llr(xi, xpi, Xc, lani,  LAnc, mc, σ²c)::Float64
        if -randexp() < llr
          llc        += llr::Float64
          Xc[xi,:]    = xpi::Array{Float64,1}
          δXc[:,:,xi] = δxi::Array{Float64,2}
          LAnc[xi,:]  = lani::Array{Float64,1}
        end
      end

    end

    return llc::Float64
  end
end



function makellr_XRupd2(δt   ::Vector{Float64}, 
                        wcol ::Array{Array{Int64,1},1})

    δt1 = δt[1]
    wci = wcol[1]

    function fx(xi  ::Int64,
        xpi ::Array{Float64,1},
        X   ::Array{Float64,2},
        lapi::Array{Float64,1},
        LA  ::Array{Float64,2},
        m   ::Float64,
        σ²  ::Float64)

        # normal likelihoods
        llr::Float64 = 0.0

        @inbounds begin
        # loop for parent nodes
            δxim1 = δt[xi-1]
            for j = wcol[xi-1]
                llr += llrdnorm_x(xpi[j], X[xi,j], 
                    X[xi-1,j] + Eδx(LA[xi-1,j], m, δxim1), 
                    δxim1*σ²)
            end

            # loop for daughter nodes
            δxi = δt[xi]
            for j = wcol[xi]
                llr += llrdnorm_μ(X[xi+1, j],
                        xpi[j]  + Eδx(lapi[j],  m, δxi),
                        X[xi,j] + Eδx(LA[xi,j], m, δxi),
                        δxi*σ²)
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
    make_mhr_upd_Xbr(Xnc1::Array{Int64,1}, Xnc2::Array{Int64,1}, wcol::Array{Array{Int64,1},1}, m::Int64, ptn::Array{Float64,1}, wXp::Array{Int64,1}, λlessthan::Int64, narea::Int64, Xupd_llr, Rupd_llr)

Make X branch DA update for a single branch using BB proposals.
"""
function make_mhr_upd_Xbr2(wcol               ::Array{Array{Int64,1},1},
                          nstep                  ::Int64,
                          ntip               ::Int64,
                          nedge              ::Int64,
                          bridx              ::Array{UnitRange{Int64},1},
                          brδt               ::Array{Array{Float64,1},1},
                          total_llr          ::Function)

  Xp   = zeros(nstep, ntip)
  δXp  = fill( NaN, ntip, ntip, nstep)
  LAnp = zeros(nstep, ntip)

  function f(br     ::Int64,
             Xc     ::Array{Float64,2},
             σ²c    ::Float64,
             mc     ::Float64,
             αc     ::Float64,
             σ²ϕ    ::Float64,
             llc    ::Float64,
             LAnc   ::Array{Float64,2},
             δXc    ::Array{Float64,3},
             δYc    ::Array{Float64,3})

    copyto!(Xp, Xc)

    upbranchX!(br, Xp, bridx, brδt, σ²ϕ)

    deltaX2!(δXp, Xp, wcol, nstep, ntip)
    sde2!(LAnp, δXp, δYc, wcol, nstep, αc, ntip)

    llr = total_llr(Xc, Xp,  LAnc, LAnp, σ²c, mc)

    if -randexp() < (llr + llr_bm(Xc, Xp, bridx[br], brδt[br], σ²ϕ))::Float64
      llc += llr::Float64
      copyto!(Xc,   Xp)
      copyto!(δXc,  δXp)
      copyto!(LAnc, LAnp)
    end

    return llc::Float64
  end
end





"""
    make_mhr_upd_Xtrio(Xnc1::Array{Int64,1}, Xnc2::Array{Int64,1}, wcol::Array{Array{Int64,1},1}, m::Int64, ptn::Array{Float64,1}, wXp::Array{Int64,1}, λlessthan::Int64, narea::Int64, Xupd_llr, Rupd_llr)

Make X trio DA update for a node and adjacent branches using BB proposals.
"""
function make_mhr_upd_Xtrio2(wcol               ::Array{Array{Int64,1},1},
                            nstep                  ::Int64,
                            ntip               ::Int64,
                            nedge              ::Int64,
                            brl                ::Array{Float64,1},
                            bridx              ::Array{UnitRange{Int64},1},
                            brδt               ::Array{Array{Float64,1},1},
                            total_llr          ::Function)

  Xp   = zeros(nstep, ntip)
  δXp  = fill( NaN, ntip, ntip, nstep)
  LAnp = zeros(nstep, ntip)

  function f(trio   ::Array{Int64,1},
             Xc     ::Array{Float64,2},
             σ²c    ::Float64,
             mc     ::Float64,
             αc     ::Float64,
             σ²ϕ    ::Float64,
             llc    ::Float64,
             LAnc   ::Array{Float64,2},
             δXc    ::Array{Float64,3},
             δYc    ::Array{Float64,3} )

    copyto!(Xp, Xc)

    pr, d1, d2 = trio

    uptrioX2!(pr, d1, d2, Xp, bridx, brδt, brl, σ²ϕ, nedge)

    deltaX2!(δXp, Xp, wcol, nstep, ntip)
    sde2!(LAnp, δXp, δYc, wcol, nstep, αc, ntip)

    
    llr = total_llr(Xc, Xp, LAnc, LAnp, σ²c, mc)

    if -randexp() < (llr + 
                     ((pr != nedge) ? 
                      llr_bm(Xc, Xp, bridx[pr], brδt[pr], σ²ϕ) : 0.0) +
                      llr_bm(Xc, Xp, bridx[d1], brδt[d1], σ²ϕ) +
                      llr_bm(Xc, Xp, bridx[d2], brδt[d2], σ²ϕ))::Float64
      llc += llr::Float64
      copyto!(Xc,     Xp)
      copyto!(δXc,   δXp)
      copyto!(LAnc, LAnp)
    end

    return llc::Float64
  end

end

"""
    Xupd_linavg!(δxi  ::Array{Float64,2},
                 lapi ::Array{Float64,1},
                 lani ::Array{Float64,1},
                 ldi  ::Array{Float64,2},
                 wci  ::Array{Int64,1},
                 xpi  ::Array{Float64,1},
                 xi   ::Int64,
                 xj   ::Int64,
                 Y    ::Array{Int64,3},
                 δyi  ::Array{Float64,2},
                 narea::Int64)

Re-estimate lineage specific means 
for a node update when `ωx >= 0.0`
"""
function Xupd_linavg2!(δxi  ::Array{Float64,2},
                       δYc  ::Array{Float64,3},
                      lani ::Array{Float64,1},
                      wcol ::Array{Array{Int64,1},1},
                      xpi  ::Array{Float64,1},
                      xi   ::Int64,
                      xj   ::Int64,
                      m    ::Float64,
                      αc    ::Float64)

  @inbounds begin
    # estimate pairwise distances
    wci = wcol[xi]
    for l = wci, j = wci
      l == j && continue
      δxi[j,l] = xpi[j] - xpi[l]
    end

    # estimate lineage averages
    for l = wci 
      lani[l] = 0.0
      for j = wci
        j == l && continue
        y = δYc[j,l,xi]
        iszero(y) && continue
        x = δxi[j,l]
        lani[l] += sign(x) * m * exp(-αc*x^2)
      end
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

function deltaX2!(δX  ::Array{Float64,3}, 
                 X    ::Array{Float64,2},
                 wcol ::Array{Array{Int64,1},1},
                 nstep::Int64,
                 ntip ::Int64)

  @inbounds begin

    for i = Base.OneTo(nstep), l = wcol[i], j = wcol[i]
      l == j && continue
      δX[j,l,i] = X[i,j] - X[i,l]
    end

  end

  return nothing
end

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
function uptrioX2!(pr   ::Int64, 
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


function mhr_upd_σ²2(σ²c      ::Float64,
                    Xc       ::Array{Float64,2},
                    llc      ::Float64,
                    prc      ::Float64,
                    σ²tn     ::Float64,
                    LAc      ::Array{Float64,2},
                    mc       ::Float64,
                    σ²prior  ::Float64,
                    σ²upd_llr::Function)

  σ²p = mulupt(σ²c, rand() < 0.3 ? σ²tn : 4.0*σ²tn)::Float64

  #likelihood ratio
  llr = σ²upd_llr(Xc, LAc, mc, σ²c, σ²p)::Float64

  # prior ratio
  prr = llrdexp_x(σ²p, σ²c, σ²prior)

  if -randexp() < (llr + prr + log(σ²p/σ²c))
    llc += llr::Float64
    prc += prr::Float64
    σ²c  = σ²p::Float64
  end

  return (llc, prc, σ²c)::Tuple{Float64,Float64,Float64}
end

function mhr_upd_α( αc      ::Float64,
                    Xc      ::Array{Float64,2},
                    δXc     ::Array{Float64,3},
                    δYc     ::Array{Float64,3},
                    llc     ::Float64,
                    prc     ::Float64,
                    αtn     ::Float64,
                    LAc     ::Array{Float64,2},
                    wcol     ::Array{Array{Int64,1},1},
                    mc      ::Float64,
                    σ²c     ::Float64,
                    αprior  ::Float64,
                    nstep   ::Int64,
                    αupd_llr::Function)
#=  αp = mulupt(αc, rand() < 0.3 ? αtn : 4.0*αtn)::Float64=#
  log_αp = StrawHat(log(αc), αtn)::Float64
  αp = exp(log_αp)        #likelihood ratio
  
  llr,LAp = αupd_llr(Xc, δXc, δYc, LAc, wcol, mc, αp, σ²c, nstep)
  # prior ratio
  prr = llrdexp_x(αp, αc, αprior)
  if -randexp() < (llr + prr + log(αp/αc))
    llc += llr::Float64
    prc += prr::Float64
    αc  = αp::Float64
    copyto!(LAc, LAp)
  end
  return (llc, prc, αc, LAc)::Tuple{Float64,Float64,Float64, Matrix{Float64}}
end




tuning_scaler(tn ::Float64, ar ::Float64) = tn*tan(ar*π/2)/tan(0.3*π/2)::Float64














function mhr_upd_m( mc      ::Float64,
                    Xc      ::Array{Float64,2},
                    llc     ::Float64,
                    prc     ::Float64,
                    mtn     ::Float64,
                    LAc     ::Array{Float64,2},
                    σ²c     ::Float64,
                    mprior  ::Float64,
                    nstep   ::Int64,
                    mupd_llr::Function)
  log_mp = StrawHat(log(mc), mtn)::Float64
        #likelihood ratio
  mp = exp(log_mp)
  llr = mupd_llr(Xc, LAc, mc, mp, σ²c)::Float64
  #if llr==0.0
  #end
  # prior ratio
  prr = llrdexp_x(mp, mc, mprior)

  if -randexp() < (llr + prr+log(mp/mc))
    llc += llr
    prc += prr
    mc  = mp
  end
  return (llc, prc, mc)::Tuple{Float64,Float64,Float64}
end







function deltaX2!(δX   ::Array{Float64,3}, 
 X    ::Array{Float64,2},
 wcol ::Array{Array{Int64,1},1},
 nstep    ::Int64,
 ntip ::Int64)

@inbounds begin

    for i = Base.OneTo(nstep), l = wcol[i], j = wcol[i]
      l == j && continue
      δX[j,l,i] = X[i,j] - X[i,l]
  end

end

return nothing
end

function sde2!(LAn    ::Array{Float64,2},
   δX     ::Array{Float64,3},                   
   δY   ::Array{Float64,3},
   wcol   ::Array{Array{Int64,1},1},
   nstep  ::Int64,
   α      ::Float64,
   ntip   ::Int64)

@inbounds begin

    for i = Base.OneTo(nstep), j = wcol[i]
      LAn[i,j] = 0.0
      for l = wcol[i]
        l == j && continue
        y = δY[l,j,i]
        iszero(y) && continue
        x = δX[l,j,i]
        LAn[i,j] += sign(x) * exp(-α*x^2)
    end
end
end

return nothing
end

function deltaXY!(δX   ::Array{Float64,3}, 
  δY   ::Array{Float64,3},
  X    ::Array{Float64,2},
  Y    ::Array{Int64,3},
  wcol ::Array{Array{Int64,1},1},
  nstep    ::Int64,
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


function makellf2(δt   ::Array{Float64,1},
  Y    ::Array{Int64,3}, 
  ntip ::Int64, 
  nstep::Int64,
  nedge::Int64)

    wf23 = Int64[]
    for j = Base.OneTo(ntip)
        push!(wf23, findfirst(Y[:,j,1] .!= 23))
    end
    wl23 = Int64[]
    for j = Base.OneTo(ntip)
        lastindex = 0
        while (wf23[j]+lastindex<nstep) &&  Y[wf23[j]+lastindex+1, j,1].!=23
            lastindex+=1
        end
        push!(wl23, wf23[j]+lastindex-1)
    end
    # number of normal evaluations
    n = 0
    for i = wf23
        n += length(i:(nstep-1))
    end
    # normal constant
    normC = -0.5*log(2.0π)
    function llf(X     ::Array{Float64,2},
        LA    ::Array{Float64,2},
        σ²    ::Float64,
        m     ::Float64)

        ll::Float64 = normC
        @inbounds begin

            # trait likelihood
            for j = Base.OneTo(ntip)
                @simd for i = wf23[j]:wl23[j]
                    ll += logdnorm_tc(X[(i+1),j], 
                        X[i,j] + Eδx(LA[i,j], m, δt[i]), 
                        δt[i]*σ²)::Float64

                end
            end
        end
        return ll::Float64
    end
    function llrf(Xc     ::Array{Float64,2},
        Xp     ::Array{Float64,2}, 
        LAc    ::Array{Float64,2},
        LAp    ::Array{Float64,2},
        σ²     ::Float64,
        m     ::Float64)

        llr::Float64 = 0.0

        @inbounds begin

        # trait likelihood
            for j = Base.OneTo(ntip)
                @simd for i = wf23[j]:(wl23[j])
                    llr += llrdnorm_xμ(Xp[(i+1),j], Xc[(i+1),j],
                        Xp[i,j] + Eδx(LAp[i,j], m, δt[i]), 
                        Xc[i,j] + Eδx(LAc[i,j], m, δt[i]),
                        δt[i]*σ²)::Float64
                end
            end
        end
        return llr::Float64
    end
    return llf, llrf
end

llrdnorm_xμ(xp::Float64, xc::Float64, μp::Float64, μc::Float64, σ²::Float64) = ((xc - μc)^2 - (xp - μp)^2)/(2.0σ²)





function makellr_σmαupd(δt  ::Vector{Float64}, 
    Y   ::Array{Int64,3}, 
    ntip::Int64)


    # which is not 23 (i.e., NaN) in each column
    w23 = UnitRange{Int64}[]
    for i = Base.OneTo(ntip)
        non23 = findall(!isequal(23), Y[:,i,1])
        push!(w23, non23[1]:non23[end-1])
    end

    function fσ(X  ::Array{Float64,2},
            LA ::Array{Float64,2},
            mc ::Float64,
            σ²c::Float64,
            σ²p::Float64)
        llr::Float64 = 0.0

        @inbounds begin
            for j = Base.OneTo(ntip)
                @simd for i = w23[j]
                    llr += llrdnorm_σ²(X[(i+1),j],  
                            X[i,j] + Eδx(LA[i,j], mc, δt[i]), 
                            δt[i]*σ²p, δt[i]*σ²c)
                    
                end
            end

        end
        return llr::Float64
    end

    function fm(X  ::Array{Float64,2},
        LA ::Array{Float64,2},
        mc ::Float64,
        mp::Float64,
        σ²c::Float64)

        llr::Float64 = 0.0

        @inbounds begin
            for j = Base.OneTo(ntip)
                @simd for i = w23[j]
                    llr += llrdnorm_μ(X[(i+1),j], 
                            X[i,j] + Eδx(LA[i,j], mp, δt[i]), #μp 
                            X[i,j] + Eδx(LA[i,j], mc, δt[i]), #μc 
                            δt[i]*σ²c)
                end
            end
        end
        return llr::Float64
    end

    function fα(X::Array{Float64,2},
        δX ::Array{Float64,3},
        δY ::Array{Float64,3},
        LAc ::Array{Float64,2},
        wcol::Array{Array{Int64,1},1},
        mc ::Float64,
        αp ::Float64,
        σ²c::Float64,
        nstep::Int64)

        llr::Float64 = 0.0
        LAp = fill(NaN, nstep, ntip)

        for i = Base.OneTo(nstep), j = wcol[i]
            LAp[i,j] = 0.0
            for l = wcol[i]
                l == j && continue
                y = δY[l,j,i]
                iszero(y) && continue
                x = δX[l,j,i]
                LAp[i,j] += sign(x) * exp(-αp*x^2)
            end
        end
        @inbounds begin
            for j = Base.OneTo(ntip)
                @simd for i = w23[j]
                    llr += llrdnorm_μ(X[(i+1),j], 
                            X[i,j] + Eδx(LAp[i,j], mc, δt[i]), #μp 
                            X[i,j] + Eδx(LAc[i,j], mc, δt[i]), #μc 
                            δt[i]*σ²c)
                end
            end
        end
        return llr::Float64,LAp::Array{Float64,2}
    end
    return fσ, fm, fα
end

##New kernel

function Bactrian(p ::Float64, tn ::Float64, m::Float64=0.95)
  s = √(1-m^2)
  z = m + randn() + s
  rdunif = rand() < 0.5
  sign = rdunif ? -1 : 1
  z=z*sign
  return z
end

function StrawHat(p::Float64, tn::Float64, a::Float64=1., b::Float64 = 1.35)
  y=.0
  u1, u2,u3 = rand(3)
  u1<a/(3*b-a) ? y = a*u2^(1/3) : y = rand(Uniform(a,b))
  if u3<1/2
   y=-y
  end 
  return p + tn*y
end


"""
    addupt(p::Float64, tn::Float64)

Gaussian parameter window move to vector.
"""
function addupt2!(p::Vector{Float64}, j::Int64, σ²::Float64)

  @inbounds p[j] += randn() * σ²

  return nothing
end
