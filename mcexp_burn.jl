#=

Burning phase for tribe.

Ignacio Quintero Mächler

t(-_-t)

April 27 2017

=#




"""
    burn_compete(...)

Burning & adaptive phase for MCMC.
"""
function burn_mcexp(total_llf     ::Function,
                    Xupd_llr      ::Function,
                    Rupd_llr      ::Function,
                    σ²upd_llr     ::Function,
                    mupd_llr     ::Function,
                    αupd_llr     ::Function,
                    mhr_upd_Xbr   ::Function,
                    mhr_upd_Xtrio ::Function,
                    nedge   ::Int64,
                    Xc      ::Array{Float64,2},
                    Yc      ::Array{Int64,3},
                    δXc     ::Array{Float64,3},
                    δYc     ::Array{Float64,3},
                    LAnc    ::Array{Float64,2},
                    σ²c     ::Float64,
                    αc      ::Float64,
                    mc      ::Float64,
                    Xnc1    ::Array{Int64,1},
                    Xnc2    ::Array{Int64,1},
                    brl     ::Array{Float64,1},
                    wcol    ::Array{Array{Int64,1},1},
                    bridx_a ::Array{Array{UnitRange{Int64},1},1},
                    trios   ::Array{Array{Int64,1},1},
                    wXp     ::Array{Int64,1},
                    σ²prior ::Float64,
                    αprior ::Float64,
                    mprior ::Float64,
                    np      ::Int64,
                    parvec  ::Array{Int64,1},
                    nburn   ::Int64,
                    screen_print::Int64,
                    obj_ar  ::Float64   = 0.234,
                    tune_int::Int64     = 1000)

  nstep, ntip, narea = size(Yc)
  # likelihood and prior
  
  llc = total_llf(Xc, LAnc, σ²c, mc)
  prc = logdexp(  σ²c, σ²prior)+
        logdexp(  mc, mprior)+
        logdexp(  αc, αprior)

  # make scaling function
  scalef = makescalef(obj_ar)

  # rest of tuning parameters
  ptn = fill(0.1, np) 

  # initialize acceptance log
  ltn = zeros(Int64, np)
  lup = zeros(Int64, np)
  lac = zeros(Int64, np)

  # row i proposals for X
  xpi  = fill(NaN, ntip)             # proposal x slice
  δxi  = fill(NaN, ntip, ntip)       # lineage pairwise differences
  lani = fill(NaN, ntip)             # lineage average for ωx < 0

  rj = wcol[1][2]                    #rj correspond à la lignée soeur de la racine.

  # progress bar
  p = Progress(nburn, dt=screen_print, desc="burn...", barlen=20, color=:green)

  # print number of parameters
  printstyled(
    "\nσ² updates per iter = ", lastindex(filter(x -> x == 1,parvec)),
    "\nm updates per iter = ", lastindex(filter(x -> x == 2,parvec)),
    "\nα updates per iter = ", lastindex(filter(x -> x == 3,parvec)),
    "\nθ  updates per iter = ", length(parvec), "\n", color = :green)

  Xcidx = CartesianIndices(Xc)

  #start burnin
  for it = Base.OneTo(nburn)
    # Update vector
    shuffle!(parvec)
    for up = parvec
      # update X
      if up > 3 

        # X updates
        @inbounds begin

          upx = wXp[up - 3] 

          # if root
          if upx == 1
            # allocate
            @simd for i = Base.OneTo(ntip)
              xpi[i] =  Xc[1,i]
            end

            # update xi
            addupt!(xpi, ptn, 1, up)  #on update tous les xpi par randn()*ptn[up].

            xpi[rj] = xpi[1]::Float64  #on met la même valeur pour les deux branches à la racine.

            llr = Rupd_llr(xpi, Xc, σ²c)::Float64

            if -randexp() < llr
              llc     += llr::Float64
              Xc[1,:]  = xpi::Array{Float64,1}
              lac[up] += 1   # log acceptance
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
            # addupt2!(xpi, xj, σ²) 
            # xpi[xj] = rand(Normal(xppi+Eδx(LAnc[xi-1,xj], m, δt[xi]), δt[xi]σ²c))
            # update xi
             addupt!(xpi, ptn, xj, up)

            if in(upx, Xnc1)        # if an internal node
              xpi[Xcidx[Xnc2[findfirst(isequal(upx),Xnc1)]][2]] = xpi[xj]
            end

            # calculate new averages
            Xupd_linavg2!(δxi, δYc,lani, wcol, 
                         xpi, xi, xj, mc, αc)

            
            llr = Xupd_llr(xi, xpi, Xc, lani, LAnc,mc, σ²c)::Float64


            if -randexp() < llr
              llc        += llr::Float64
              Xc[xi,:]    = xpi::Array{Float64,1}
              δXc[:,:,xi] = δxi::Array{Float64,2}
              LAnc[xi,:]  = lani::Array{Float64,1}
              lac[up]    += 1   # log acceptance
            end
          end
        end
      end
      # if σ² is updated
      if up == 1

        σ²p = mulupt(σ²c, ptn[1])::Float64
        #likelihood ratio
        llr = σ²upd_llr(Xc, LAnc, mc, σ²c, σ²p)::Float64
        # prior ratio
        prr = llrdexp_x(σ²p, σ²c, σ²prior)

        if -randexp() < (llr + prr + log(σ²p/σ²c))
          llc += llr::Float64
          prc += prr::Float64
          σ²c  = σ²p::Float64
          lac[1] += 1
        end      
      #update α
      #=elseif up == 2
         log_αp = StrawHat(log(αc), ptn[2])::Float64
         αp = exp(log_αp)
      #   #likelihood ratio
         llr,LAp = αupd_llr(Xc, δXc, δYc, LAnc, wcol, mc, αp, σ²c, nstep)

      #   # prior ratio
         prr = llrdexp_x(αp, αc, αprior)

         if -randexp() < (llr + prr+ log(αp/αc))
           llc += llr::Float64
           prc += prr::Float64
           αc  = αp::Float64
           copyto!(LAnc, LAp)
           lac[2] += 1
         end=#
       # update m
#=      elseif up == 3
         log_mp = StrawHat(log(mc), ptn[3])::Float64
         mp = exp(log_mp)
         #likelihood ratio
         llr = mupd_llr(Xc, LAnc, mc, mp, σ²c)::Float64

         # prior ratio
         prr = llrdexp_x(mp, mc, mprior)

         if -randexp() < (llr + prr+ log(mp/mc))
           llc += llr
           prc += prr
           mc  = mp
           lac[3] += 1 
         end=#
      end

      # log number of updates
      ltn[up] += 1
      lup[up] += 1
      if (in(tune_int,ltn))
        wts = findall(isequal(tune_int),ltn)      # which to scale
        for j = wts
          ar     = lac[j]/lup[j]
          ptn[j] = tuning_scaler(ptn[j],ar)
          ltn[j] = 0
          lac[j] = 0
          lup[j] = 0
        end
      end

    end

    next!(p)
  end
  return llc, prc, Xc, Yc, LAnc, δXc, δYc, σ²c, αc, mc, ptn
end
