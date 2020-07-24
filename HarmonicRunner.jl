using DelimitedFiles
using Random
using  PyCall
using LinearAlgebra
using SparseArrays
using Pkg; Pkg.add("PyPlot")
import Pkg; Pkg.add("GraphRecipes")
using GraphRecipes
using Plots
nx=pyimport("networkx")

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
    println(typeof(cliq))
    clique = cliq[:,:]
    println(typeof(clique))
    clique = hcat(clique...)
    n = maximum(clique)[1]
    A = sparse(clique[1,:],
              clique[2,:] ,1, n, n)
    for i=1:size(A,1)
        for j=1:size(A,1)
            if A[i,j] != 0
                A[i,j] = 1
            end
        end
    end
    return A
end

h=getA()

function getA()
    hypergraph = convertToHype("email-Enron/email-Enron-nverts.txt", "email-Enron/email-Enron-simplices.txt")
    clique = hype2Clique(hypergraph)
    A1 = cliqueToAdjMatrix(clique)
    ### Calculate Centrality Measure by converting clique adjacency matrix
    ### into graph
    A = Matrix(A1)
    G = nx.Graph(A)
    return hypergraph
end
#----------------------------------------------------------------#
######## ERROR: AssertionError: length(node_weights) == n ########
#----------------------------------------------------------------#
function graphNetwork(network::SparseMatrixCSC{Int64, Int64})
    #=
    Returns: graph of clique with node size scaled by
            their harmonic centrality score
    =#
    Random.seed!(123)  # plot same each time
    G = nx.Graph(A)
    harmCent=nx.harmonic_centrality(G)
    cent=[harmCent[i] for i in 0:147]
    graphplot(network,
              markersize = 0.2,
              markercolor = "red",
              fontsize = 10,
              linecolor = :darkgrey,
              node_weights = cent.+1
              )
    savefig("one.png")
end
graphNetwork(getA())


for i in 1:length(h)
    if length(findall(h[i],1))>1
        println(i)
    end
end
