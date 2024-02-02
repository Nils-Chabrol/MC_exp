

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
        end <-begin <- rep(0, nrow(.subset2(tree,'edge'))) 

# Do it recursively

fx <- function(tree, node, begin, end){
  
  cur.time <- 0
  root <- length(.subset2(tree,'tip.label')) + 1
  if (node > root){
    cur.time <- end[which(.subset2(tree,'edge')[,2] == node)] 
  }
  
  dset <- .subset2(tree,'edge')[,2][.subset2(tree,'edge')[,1] == node]
  
  i1 <- which(.subset2(tree,'edge')[,2] == dset[1])
  i2 <- which(.subset2(tree,'edge')[,2] == dset[2])
  
  end[i1] <- cur.time + .subset2(tree,'edge.length')[i1]
  end[i2] <- cur.time + .subset2(tree,'edge.length')[i2]
  
  if (dset[1] > length(.subset2(tree,'tip.label'))){
    begin[.subset2(tree,'edge')[,1] == dset[1]] <- end[i1]
    data <- fx(tree, node = dset[1], begin, end)
    tree<-data[[1]]
    begin<-data[[2]]
    end<-data[[3]]  
    
  }
  if (dset[2] > length(.subset2(tree,'tip.label'))){
    begin[.subset2(tree,'edge')[,1] == dset[2]] <- end[i2]
    data <- fx(tree, node = dset[2], begin, end)
    tree<-data[[1]]
    begin<-data[[2]]
    end<-data[[3]]
    
  }
  
  return(list(tree, begin, end))
  
}

data <- fx(tree, node = length(.subset2(tree,'tip.label')) + 1, begin, end)
tree<-data[[1]]
begin<-data[[2]]
end<-data[[3]]

maxbt <- max(end)
nodes <- (length(.subset2(tree,'tip.label')) + 1):(2*length(.subset2(tree,'tip.label')) - 1)
bt <- numeric(length(nodes))
names(bt) <- nodes
for (i in 1:length(bt)){
  tt <- begin[.subset2(tree,'edge')[,1] == nodes[i]][1]
  bt[i] <- maxbt - tt
}
bt


        """)
    brtimes = rcopy(brtimes)
    return tree, brtimes
else
    return tree
end
end


function read_data_mcexp(tree_file::String,
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

  tip_values = Dict(tip_labels[val] => data_values[i] 
  for (i,val) = enumerate(data_tlab)) #attribute a value to a tip number
    tip_areas = Dict(tip_labels[val] => data_areas[i,:] 
    for (i,val) = enumerate(data_tlab))
      return tip_values, tip_areas, tree, bts
end




"""
    initialize_data(tip_values::Dict{Int64,Float64}, tip_areas ::Dict{Int64,Array{Int64,1}}, m::Int64, tree::rtree, bts::Array{Float64,1})

Function to initialize `X` and `Y` matrix given
the tip_values and tip_areas (as Dictionaries).
"""
function initialize_data_mcexp(tip_values::Dict{Int64,Float64},
 tip_areas ::Dict{Int64,Array{Int64,1}},
 min_dt    ::Float64,
 tree      ::rtree,
 bts       ::Array{Float64,1})

  br     = branching_times(tree) # parent node ; child node ; 3 : branch length ; 4 : parent height ; 5 : children height.
  n      = tree.nnod + 1
  nareas = length(tip_areas[1])
  bts = unique(br[:,4])
  #*
  # make times and δt vector
  #*

    # make sure each branch has nareas + 1 sampling points
  ets = unique(cat(br[:,4],br[:,5], dims=(1)))
  #stocke les temps des époques #créer des nouvelles époques en fonction de la biogéographie.

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

  #create δt vector
  δt = abs.(diff(ets))

  # initialize data augmentation matrices
  X = fill(NaN, length(ets), n) #each lign correspond to the trait value at one point in time, each column represents a point in time.
  B = copy(X)
  Y = fill(23, length(ets), n, nareas)

  # which rows are branching points
  wch = indexin(bts, ets)

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
  

  Extinct=Int64[]
  for i in eachindex(br[:,1])
    if br[i,2]<=n && br[i,5]!=0.0
      push!(Extinct, i)
    end
  end

  Extinction_times = Array{Real}(undef,length(Extinct), 3) #Extinction_times : [branch that goes extinct ; time of extinction ; position of the time of extinction in ets.]
  for i in 1:lastindex(Extinct)
    Extinction_times[i,1] = Int64(br[Extinct[i],2])
    Extinction_times[i,2] = br[Extinct[i],5]
    Extinction_times[i,3] = indexin(Extinction_times[i,2], ets)[1]
  end


  ord=sortperm(br[:,5])
  ord2 = filter(x -> !(br[x,5] in Extinction_times[:,2]), ord)
  bord = ord2[(1:tree.nnod-1) .+ n.-length(Extinct)]

  for i = Base.OneTo(tree.nnod-1) #retrouve les noeuds créer à t0
    fn  = br[bord[i],2] #first node
    wda = br[findall(isequal(fn), br[:,1]),2]

    cda = indexin(wda,alive)
    coup[i,1:2] = cda

    alive[cda] = [fn,0]
    wrtf = wts[i+1]:(wts[i]-1)

    B[wrtf,:] = repeat(alive, inner = length(wrtf))

    wca = findall(x -> x > 0, alive)

    X[wts[i+1]:length(ets),wca] .= 1.0
  end

  for i in 1:length(Extinct)
    tax_ind = Extinction_times[i,1]
    time_ind = (Extinction_times[i,3]+1):length(ets)
    B[time_ind,tax_ind].=X[time_ind,tax_ind].=NaN
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
  si = initialize_X!(tip_values, X, B, Extinction_times, ncoup, δt, tree)
  initialize_Y!(tip_areas, Y, B, Extinction_times)

  # declare non-23s for Y

  return X, Y, B, ncoup, δt, tree, si
end


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
    @views tree_depth = maximum(brs[:,5])

    for j = eachindex(tree.el) 
      brs[j,4] = tree_depth - brs[j,4]
      brs[j,5] = tree_depth - brs[j,5]
      if brs[j,5]<.00001
        brs[j,5] = 0.0
      end
    end
  end
  return brs
end


function initialize_X!(tip_values::Dict{Int64,Float64},
                       X         ::Array{Float64,2},
                       B         ::Array{Float64,2},
                       Extinction_times ::Array{Real,2},
                       ncoup     ::Array{Int64,2},
                       δt        ::Array{Float64,1},
                       tree      ::rtree)

  co1    = ncoup[:,1]
  co2    = ncoup[:,2]
  nr, nc = size(X)
  nbrs   = 2*nc - 1 #number of branches

  # matrix of delta times
  δtM = zeros(lastindex(δt)+1,nc)
  for i = Base.OneTo(nc)
    δtM[2:end,i] = δt
  end   #δtM représente une matrice avec tous les pas de temps.

  Xnod = ncoup[size(ncoup,1):-1:1,:] #renverse l'ordre des lignes de ncoup.

  sort!(Xnod, dims=1) #trie Xnod des valeurs les plus basses (ligne 1) aux valeurs les plus hautes.

  # brownian motion
  bm_ll = make_bm_ll(tip_values, tree) #construit un fonction de ll pour le bm
  op    = optimize(bm_ll, rand(nc), Options(g_tol = 1e-6, iterations = 100_000, store_trace = false)) #on optimise les valeurs aux noeuds selon un bm sur l'arbre en fonction des valeurs aux tips.

  ar = minimizer(op)[2:end] 
  si = minimizer(op)[1] #valeur optimal de σ

  X[Xnod[:,1]] = ar #on place les valeurs aux noeuds dans X.

  for i in Base.OneTo(nc)
    if i in Extinction_times[:,1]
      ind = indexin(i, Extinction_times[:,1])[1]
      X[Extinction_times[ind,3],i] = tip_values[i]
    else
      X[nr,i] = tip_values[i]
    end #on place les valeurs aux tips dans X.
  end

  X[co2] = X[co1] #les noeuds appairés ont les mêmes valeurs.

  for i = setdiff(1:nbrs,nc+1) #pour tous les numéros de branches sauf la racine

    wbranch = findall(isequal(i), B) #Dans B sont stockés les numéros des branches, l'arbre est représenté sous forme de matrice qui permet de bien indexé les éléments de X.
    l_i = wbranch[1] - CartesianIndex(1,0)
    l_f = wbranch[end]
    X[CartesianIndices((l_i[1]:l_f[1], l_i[2]:l_i[2]))] = 
      bb(X[l_i], X[l_f], δtM[wbranch])
  end

  X[co2] .= NaN #on enlève les valeurs avant la spéciation.

  return si
end

function initialize_Y!(tip_areas::Dict{Int64,Array{Int64,1}},
                       Y        ::Array{Int64,3},
                       B        ::Array{Float64,2},
                       Extinction_times ::Array{Real,2})

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
    if i in Extinction_times[:,1]
      indice = indexin(i, Extinction_times[:,1])[1]

      Y[Extinction_times[indice,3],i,:] = tip_areas[i]
    else
      Y[end,i,:] = tip_areas[i]
    end
  end
end