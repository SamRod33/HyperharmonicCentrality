using DelimitedFiles
using Random
using  PyCall
using LinearAlgebra
using SparseArrays
using Pkg; Pkg.add("PyPlot")
import Pkg; Pkg.add("GraphPlot")
using GraphPlot
using Plots
using LightGraphs
using Colors
using Cairo
using Compose
using PyPlot
nx=pyimport("networkx")


function get_val(dict::Dict{Int64,Int64}, key::Int64)
    if !haskey(dict, key)
        n = length(dict) + 1
        dict[key] = n
        return n
    end
    return dict[key]
end



function convertToHype(dataFile1::String, dataFile2::String)
    #=
    Returns: clique representing the projected graph of
                a hypergraph

    Parameter: File containing hypergraph
    Precondition: files exist and are in the same folder, first file lists the
    hyperedges and second file lists the nodes
    =#
    #FileReading Hypergraph (if we have two Files)

    #first file's line number* is its position in the matrix
    #the number on that line* is how many nodes are in the edge
    #go into second file and take that many nodes from last (keep track of
    #position last taken from)

    #final array like this:
    # hyper = [ [nodes in hyperedge1]; [nodes in hyperedge2]; ... ;[nodes in hyperedge n] ]

    edges = readdlm(dataFile1)
    nodes = readdlm(dataFile2)
    numedges = size(edges,1) #number of hyperedges / rows in hypergraph matrix
    numnodes = size(nodes,1) #total number of hypernodes

    #initialize hypergraph as an array of arrays
    hyper = hyper=[[1,2]]

    #go through each line
    i = 1 #line pos we're at in second file
    for line in 1:numedges
        take = edges[line] #number of nodes to take from other file
        if take > 1
            Edge = nodes[i:(i-1)+Int64(take)] #makes an array representing hyperedge with its nodes
            push!(hyper, Edge) #place hyperedge array in overall hypergraph
        end
        i += Int64(take) #update line postion in second file
    end
    hyper=hyper[2:end]

    # Hypergraph (Converting if hyperedge is line by line)
    # line2vec(line) = parse.(Int, split(line, '\t' ))
    # hyper=[line2vec(line) for line in eachline(dataFile;
    #         keep = false) if !isempty(line)]

    node_map = Dict{Int64, Int64}();
    for i in 1:length(hyper)
        for j in 1:length(hyper[i])
            hyper[i][j]=get_val(node_map,hyper[i][j])
        end
    end

    return unique(hyper)
    # Once we have clique, this makes the adjency matrix of clique
 end


function hype2Clique(hyper::Vector{Array{Int64, 1}})
    # Conversion to clique
    clique = Array[]

    m = size(hyper, 1) # amount of hyperedges
    for hyperedge = 1:m
        n = size(hyper[hyperedge], 1) # nodes in the current hyperedge
        for posNodeI=1:n
            for posNodeJ=posNodeI+1:n
                i = hyper[hyperedge][posNodeI] # current node
                j= hyper[hyperedge][posNodeJ] # running through all
                                                 # other nodes in hyperedge
                push!(clique, [i,j]) # create pair (i,j)
            end
        end
    end

    return unique(clique)
end


function cliqueToAdjMatrix(cliq::Vector{Array})
    # Converts Clique to Adjacency Matrix representation
    # Makes all values of weights to be one
    clique = cliq[:,:]
    clique = hcat(clique...)
    n = maximum(clique)[1]
    A = sparse(clique[1,:],
              clique[2,:] ,1, n, n)
    A = A+ A'
    for i=1:size(A,1)
        for j=1:size(A,1)
            if A[i,j] != 0
                A[i,j] = 1
            end
        end
    end
    return A
end


function cliqueToAdjMatrix_Weighted(cliq::Vector{Array})
    # Converts Clique to Adjacency Matrix representation
    # Makes all values of weights to be one
    clique = cliq[:,:]
    clique = hcat(clique...)
    n = maximum(clique)[1]
    A = sparse(clique[1,:],
              clique[2,:] ,1, n, n)
    A = A+ A'
    return A
end


function getHypergraph()
    #= Returns Hypergraph from files Enron =#
    hypergraph = convertToHype("email-Enron/email-Enron-nverts.txt", "email-Enron/email-Enron-simplices.txt")
    return hypergraph
end


function getClique()
    #= Returns clique from hypergraph of Enron=#
    hypergraph = getHypergraph()
    clique = hype2Clique(hypergraph)
    return clique
end


function getAdjMatrix()
    clique = getClique()
    A1 = cliqueToAdjMatrix(clique)
    ### Calculate Centrality Measure by converting clique adjacency matrix
    ### into graph
    A = Matrix(A1)
    return A
end


function getAdjMatrix_Weighted()
    #= Returns clique that is supposed to be weighted =#
    hypergraph = getHypergraph()
    clique = hype2Clique(hypergraph)
    A1 = cliqueToAdjMatrix_Weighted(clique)
    A = Matrix(A1)
    return A
end


function getHarmonicCentrality(network::Array{Int64, 2})
    #= Returns harmonic centrality score for each node as a dictionary =#
    G = nx.Graph(network)
    harmCent=nx.harmonic_centrality(G)
    return harmCent
end


function getClosenessCentrality(network::Array{Int64,2})
    #= Returns closeness centrality score for each node as a dictionary =#
    G = nx.Graph(network)
    closeCent=nx.closeness_centrality(G)
    return closeCent
end


#######################################################################
#------------------DATA for Enron Data Set:---------------------------#
#######################################################################
# nodes, edges, hyperedges
# nodes      - 143
# edges      - 2168
# hyperedges - 1487


#######################################################################
#---------------------------------------------------------------------#
#---------------------STATISTICS ANALYSIS-----------------------------#
#---------------------------------------------------------------------#
#######################################################################

#######################################################################
#--------Calculating L_2 Norm (How close are the two datasets)--------#
#######################################################################
hCent = getHarmonicCentrality(getAdjMatrix())
ordered_hCent=sort(hCent)
hCentScores = [t[2] for t in ordered_hCent]
hCentScoresNorm = hCentScores / maximum(hCentScores)

cCent = getClosenessCentrality(getAdjMatrix())
ordered_cCent=sort(cCent)
cCentScores = [t[2] for t in ordered_cCent]
cCentScoresNorm = cCentScores/ maximum(cCentScores)

norm(hCentScoresNorm-cCentScoresNorm,2)
### L_2 Norm --> 0.1869 result on a normalized scale (i.e from 0-1 scale) ###

########################################################################
#---------------------HAS YET TO WORK----------------------------------#
########################################################################
function rank_corrbetweenness(name::Any,dataset::Any, index_start::Int64)
    cec_c=dataset[:,1]
    n=length(cec_c)
    cec_c[cec_c]=[i for i in 1:n]
    hec_c=dataset[:,2]
    hec_c[hec_c]=[i for i in 1:n]
    start = log10(index_start - 1)
    finish = log10(length(cec_c) - 1)
    ran = [convert(Int64, round(v)) for v in 10 .^ range(start, stop=finish, length=500)]
    sp = sortperm(cec_c)
    cec_c, hec_c= cec_c[sp], hec_c[sp]
    cec_hec = [corspearman(cec_c[(end-s):end], hec_c[(end-s):end]) for s in ran]
    sp = sortperm(hec_c)
    cec_c, hec_c = cec_c[sp], hec_c[sp]
    hec_cec = [corspearman(hec_c[(end-s):end], cec_c[(end-s):end]) for s in ran]
    clf()
    xs = collect(ran) .+ 1
    semilogx(xs, cec_hec, linestyle="-",  lw=1,    label="TWDB-H2NB")
    semilogx(xs, hec_cec, linestyle=":",  lw=2.25, label="H2NB-TWDB")
    fsz = 20
    ax = gca()
    ax.tick_params(axis="both", labelsize=fsz-2, length=7, width=1.5)
    ax.tick_params(axis="x",    which="minor", length=4, width=1)
    legend(fontsize=fsz-10, frameon=false)
    title("$name", fontsize=fsz)
    xlabel("Number of top-ranked elements", fontsize=fsz)
    ylabel("Spearman’s rank corr. coeff.", fontsize=fsz)
    savefig("6241-spearman-$name-rankcorr.png", bbox_inches="tight")
end
rank_corrbetweenness("B", h, 10)

##____IDEA_____##
# Rank top ten closness and harmonic (normalized)
# Generalize for each rank afterwards (spearmans correlation coefficient)
##_____________##

##############################################################################
######-------------- WHAT IF WE DIDNT MAKE ALL WEIGHTS ONE ------------#######
##############################################################################
hCentW=getHarmonicCentrality(getAdjMatrix_Weighted())
ordered_hCentW=sort(hCentW)
hCentWScores = [t[2] for t in ordered_hCentW]
hCentWScoresNorm = hCentWScores / maximum(hCentWScores)

hCentNotW=getHarmonicCentrality(getAdjMatrix())
ordered_hCentNotW=sort(hCentNotW)
hCentNotWScores = [t[2] for t in ordered_hCentNotW]
hCentNotWScoresNorm = hCentNotWScores / maximum(hCentNotWScores)

norm(hCentNotWScoresNorm-hCentWScoresNorm,2)
# L_2 Norm turned out to be 0... I don't know what to make of this
#           the Adjacency Matrix for weight included had a max weight of 2
#           on some nodes, not sure if the renaming of nodes affected this...


#############################################################################
#-------------------- GRAPHING NETWORKS ------------------------------------#
#############################################################################
function graphHarmonicNetwork(network::Array{Int64, 2})
    #=
    Returns: graph of clique with node color scaled by
            their harmonic centrality score
    =#
    Random.seed!(12345678)  # plot same each time
    harmCent = getHarmonicCentrality(network)
    n = length(keys(harmCent))-1
    cent=[harmCent[i] for i in 0:n]
    cent = cent / maximum(cent)
    g = SimpleGraph(network)
    nodefill = [RGBA(i,0.0,0.0,i) for i in cent]
    layout=(args...)->spring_layout(args...; C=35)
    p= gplot(g, layout=layout, nodefillc=nodefill,
            edgestrokec=colorant"grey", nodelabel = 1:143,
            NODELABELSIZE=3, NODESIZE=0.040)
    draw(PNG("CliqueHarmonic_Weighted.png", 16cm, 16cm), p)
end


function graphClosenessNetwork(network::Array{Int64, 2})
    #=
    Returns: graph of clique with node color scaled by
            their closeness centrality score
    =#
    Random.seed!(12345678)  # plot same each time
    closeCent = getClosenessCentrality(network)
    n = length(keys(closeCent))-1
    cent=[closeCent[i] for i in 0:n]
    cent = cent / maximum(cent)
    g = SimpleGraph(network)
    nodefill = [RGBA(i,0.0,0.0,i) for i in cent]
    layout=(args...)->spring_layout(args...; C=35)
    p= gplot(g, layout=layout, nodefillc=nodefill,
            edgestrokec=colorant"grey", nodelabel = 1:143,
            NODELABELSIZE=3, NODESIZE=0.040)
    draw(PNG("CliqueClosnessRED.png", 16cm, 16cm), p)
end

graphHarmonicNetwork(getAdjMatrix())
graphClosenessNetwork(getAdjMatrix())
graphHarmonicNetwork(getAdjMatrix_Weighted())
