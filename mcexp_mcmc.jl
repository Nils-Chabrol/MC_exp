
function mcexp_mcmc(Xc          ::Array{Float64,2},
    Yc          ::Array{Int64,3},
    ncoup       ::Array{Int64,2},
    δt          ::Array{Float64,1},
    edges       ::Array{Int64,2},
    brl         ::Array{Float64,1},
    B           ::Array{Float64,2};
    niter       ::Int64             = 500_000,
    nthin       ::Int64             = 1_000,
    nburn       ::Int64             = 500_000,
    saveXY      ::Tuple{Bool,Int64} = (false, 1_000),
    saveDM      ::Tuple{Bool,Int64} = (false, 1_000),
    σ²prior     ::Float64           = 1e-1,
    αprior      ::Float64           = 1e-1,
    mprior      ::Float64           = 1e-1,
    out_file    ::String            = "tribe_results",
    weight      ::NTuple{4,Float64} = (0.3,0.3,0.3,5e-3),
    σ²i         ::Float64           = 1.,
    αi          ::Float64           = .5,
    mi          ::Float64           = 1.,
    stbrl       ::Float64           = 1.,
    screen_print::Int64             = 5)

    printstyled("Data successfully processed", bold = true,color=:green)
    # dims
    nstep, ntip, narea  = size(Yc)

    # coupled nodes for X
    Xnc1 = ncoup[:,1]
    Xnc2 = ncoup[:,2]

    # which nodes are not NaN in Xc, without the final tips.
    wXp = setdiff(LinearIndices(Xc)[findall(!isnan, Xc)], nstep:nstep:length(Xc))

    # tie trait coupled nodes
    Xc[Xnc2] = Xc[Xnc1]

    #create object with column indices
    wcol = create_wcol(Xc)

    # add a branch as long as the tree
    stbrl = isone(stbrl) ? sum(δt) : stbrl

    # add long stem branch 
    edges = cat(edges, [2*ntip ntip + 1], dims = 1) #ajoute un branche à la base de la racine

    push!(brl, stbrl)

    # number of edges
    nedge = size(edges,1) 

    # make edge triads
    trios = maketriads(edges) #fait les triades en fonction de la position de la branche dans la matrice edge.

    # make ragged array with index for each edge in Yc
    bridx = make_edgeind(edges[:,2], B, ntip)

    # expand bridx for each area
    bridx_a = Array{UnitRange{Int64},1}[]
    push!(bridx_a, bridx)

    for j = 2:narea
        bridxj = copy(bridx)
        for i in Base.OneTo(nedge)
        bridxj[i] = (bridxj[i]) .+ (j-1)*(nstep*ntip)
        end
    push!(bridx_a, bridxj)
    end

    stemevc = [[rand()] for i in 1:narea]
    brct = make_edgeδt(bridx, δt, nstep)

    ## allocate averages for X and Y
    # X and Y distance matrix
    δXc = fill(NaN, ntip, ntip, nstep)
    δYc = fill(NaN, ntip, ntip, nstep)
    #deltaXY!(δXc, δYc, Xc, Yc, wcol, nstep, ntip, narea)
    deltaX2!(δXc,Xc, wcol, nstep, ntip)

    # lineage averages
    LAnc = fill(NaN, nstep, ntip)

    sde2!(LAnc, δXc, δYc, wcol, nstep, αi, ntip)



    ## make likelihood and prior functions
    total_llf, total_llr = makellf2(δt, Yc, ntip, nstep, nedge)

    Xupd_llr,  Rupd_llr  = makellr_XRupd2(δt, wcol)
    σ²upd_llr, mupd_llr, αupd_llr = makellr_σmαupd(δt, Yc, ntip)
    # number of xnodes + α + m + σ² 
    np = length(wXp) + 3  

    # parameter update vector
    parvec = collect(1:np)  # parameter vector
    append!(parvec, fill(1,ceil(Int64,np*weight[1]))) # weight of σ²
    append!(parvec, fill(2,ceil(Int64,np*weight[2]))) # weight of α
    append!(parvec, fill(3,ceil(Int64,np*weight[3]))) # weight of m
    append!(parvec, repeat(4:np, inner = ceil(Int64,np*weight[4]))) #weight of the X.



    # create update functions for Xbr, Y, Ybr and XYbr
    mhr_upd_Xbr     = make_mhr_upd_Xbr2(wcol, nstep, ntip, nedge, 
    bridx, brct, total_llr)
    mhr_upd_Xtrio   = make_mhr_upd_Xtrio2(wcol, nstep, ntip, nedge,
    brl, bridx, brct, total_llr)
    # burning phase
    llc, prc, Xc, Yc, LAnc, δXc, δYc, σ²c, αc, mc,  ptn = burn_mcexp(total_llf,Xupd_llr, Rupd_llr, σ²upd_llr, mupd_llr, αupd_llr,
                                                                    mhr_upd_Xbr, mhr_upd_Xtrio, nedge,
                                                                    Xc, Yc, δXc, δYc, LAnc, σ²i, αi, mi,
                                                                    Xnc1, Xnc2, brl, wcol, bridx_a, 
                                                                    trios, wXp, σ²prior,αprior, mprior, np, parvec, 
                                                                    nburn, screen_print)

     # log probability of collision
    max_δt = maximum(δt)::Float64
    # log for nthin
    lthin = 0
    # variables to save X and Y 
    if saveXY[1]
        XYsav  = 0
        XYlit  = 0
        xylogs = fld(niter,saveXY[2])
        Xlog   = zeros(Float64, nstep, ntip, xylogs)
        Ylog   = zeros(Int64,   nstep, ntip, narea, xylogs)

        if saveDM[1]
        DMsav  = 0
        DMlit  = 0
        dmlogs = fld(niter,saveDM[2])
        LAlog = zeros(Float64, nstep, ntip, dmlogs)
        LDlog  = zeros(Float64, nstep, ntip, narea, dmlogs)
        end
    end

    # progress bar
    p = Progress(niter, dt=screen_print, desc="mcmc...", barlen=20, color=:green)

  # create X parameter update function
    mhr_upd_X = make_mhr_upd_X2(Xc, Xnc1, Xnc2, wcol, ptn, wXp, nstep,
                ntip, Xupd_llr, Rupd_llr)

    # write to file
    open(out_file*".log", "w")

    open(out_file*".log", "a") do f

        print(f, "iteration\tlikelihood\tprior\tσ²\t m \t α \n")

        #=
        start MCMC
        =#

        for it = Base.OneTo(niter)
            # update vector order
            shuffle!(parvec)

            for up = parvec

                # update X[i]
                if up > 3
                    llc = mhr_upd_X(up, Xc, δXc,δYc, σ²c ,mc , αc, llc, LAnc)
                
                #if σ² is updated
                elseif up == 1
                    llc, prc, σ²c = mhr_upd_σ²2(σ²c, Xc, llc, prc, ptn[1], 
                                                LAnc, mc, σ²prior, σ²upd_llr)
                #update α
#=                elseif up == 2
                   llc, prc, αc, LAnc = mhr_upd_α(αc,Xc, δXc, δYc, llc, prc, ptn[2], LAnc, wcol, mc, σ²c, αprior, nstep, αupd_llr)=#

                # update m
#=                else up==3
                   llc, prc, mc = mhr_upd_m(mc, Xc, llc, prc, ptn[2], LAnc, σ²c, mprior, nstep, mupd_llr)=#
                end
            end

            # make DA updates with some Pr
            #make X branch update
            if rand() < 1.0
               σ²ϕ = σ²ϕprop()
               llc = mhr_upd_Xbr(rand(Base.OneTo(nedge-1)),
                   Xc,  σ²c, mc, αc, σ²ϕ, llc, 
                   LAnc, δXc, δYc)
            end

            #make X trio update
            if rand() < 1.0
                σ²ϕ = σ²ϕprop()
                llc = mhr_upd_Xtrio(rand(trios),
                    Xc, σ²c, mc, αc, σ²ϕ, llc, 
                    LAnc, δXc, δYc)
            end # end parameter loop

            # log parameters
            lthin += 1
            if lthin == nthin
                print(f, it,"\t", llc, "\t", prc,"\t",
                σ²c,"\t",mc,"\t",αc,"\n")
                lthin = 0
            end

            # log X & Y
            if saveXY[1]
                XYsav += 1
                if XYsav == saveXY[2]
                    @inbounds begin
                        XYlit += 1
                        Xlog[:,:,  XYlit] = Xc
                        Ylog[:,:,:,XYlit] = Yc
                    end
                XYsav = 0
                end
            end

            # log deterministic matrices
            if saveDM[1]
                DMsav += 1
                if DMsav == saveDM[2]
                    @inbounds begin
                        DMlit += 1
                        LAlog[:,:,  DMlit] = ωxc > 0 ? LApc * ωxc : LAnc * ωxc
                        LDlog[:,:,:,DMlit] = LDc 
                    end
                DMsav = 0
                end
            end

        next!(p) 
        end# end MCMC
    end # end print loop


    if saveXY[1]
    # save X and Y as R objects
        @rput Xlog
        @rput Ylog
        @rput δt
        @rput B

    if saveDM[1]
        @rput LAlog
        @rput LDlog
        reval("""
            delta.t <- δt
            save(Xlog, Ylog, B, delta.t, LAlog, LDlog, 
            file = '$out_file.rda')
            """)
    else
        reval("""
            delta.t <- δt
            save(Xlog, Ylog, B, delta.t, file = '$out_file.rda')
            """)
    end

end

return llc, prc, σ²c, mc, αc, Xc, Yc

end
