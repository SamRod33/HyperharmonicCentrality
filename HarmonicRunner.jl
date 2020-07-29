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


function getAdjMatrix()
    clique = getClique()
    A1 = cliqueToAdjMatrix(clique)
    ### Calculate Centrality Measure by converting clique adjacency matrix
    ### into graph
    A = Matrix(A1)
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
    nodefill = [RGBA(0.0,i,0.0,i) for i in cent]
    layout=(args...)->spring_layout(args...; C=35)
    p= gplot(g, layout=layout, nodefillc=nodefill,
            edgestrokec=colorant"grey", nodelabel = 1:143,
            NODELABELSIZE=3, NODESIZE=0.040)
    draw(PNG("CliqueHarmonic.png", 16cm, 16cm), p)
end
graphHarmonicNetwork(getAdjMatrix())

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
    nodefill = [RGBA(0.0,i,0.0,i) for i in cent]
    layout=(args...)->spring_layout(args...; C=35)
    p= gplot(g, layout=layout, nodefillc=nodefill,
            edgestrokec=colorant"grey", nodelabel = 1:143,
            NODELABELSIZE=3, NODESIZE=0.040)
    draw(PNG("CliqueClosnessNorm.png", 16cm, 16cm), p)
end
graphClosenessNetwork(getAdjMatrix())

# STATISTICS ANALYSIS SUGGESTIONS
# sum of (A.-B) do line 193 instead
# L_2 Norm
# norm function

# Rank top ten closness and harmonic (normalized)
# Generalize for each rank afterwards (spearmans correlation coefficient)

# WHAT IF WE DIDNT MAKE ALL WEIGHTS ONE

norm(hcent-ccent,2)

# DATA:
# nodes, edges, hyperedges
# nodes      - 143
# edges      - 2168
# hyperedges - 1487
