#=
Main Module for the research project

This module provides a variety of methods to manipulate a hypergraph or graph
network given source files. We use the Enron Email Dataset to apply our
research. The calculation of Harmonic Centrality and Closeness Centrality are
currently supported, as that is our focus for the project.

This module transforms the dataset into an adjacency matrix that represents
the weighted and unweighted clique formed from its hypergraph.

Unweighted and Weighted representations of the dataset are supported.

The GraphVisuals.jl and StatisticalAnalysis.jl modules derive their calculations
from this main module

Author: Samuel Rodriguez   (sar325@cornell.edu)
        Monique Rampersaud (mar452@cornell.edu)
Date:   08.01.20
=#
using DelimitedFiles
using Random
using  PyCall
using LinearAlgebra
using LightGraphs
using SparseArrays
using StatsBase
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


 function convertToHype_Weighted(dataFile1::String, dataFile2::String)
     #=
     Returns: clique representing the projected graph of
                 a hypergraph with weights included

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

     # No longer want to return unique hypergraph
     return hyper
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


function hype2Clique_Weighted(hyper::Vector{Array{Int64, 1}})
    # Conversion to clique with edge weights included
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

    return clique
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
    hypergraph = convertToHype("email-Enron/email-Enron-nverts.txt",
    "email-Enron/email-Enron-simplices.txt")
    return hypergraph
end


function getHypergraph_Weighted()
    #= Returns Hypergraph from files Enron Weighted Version =#
    hypergraph = convertToHype_Weighted("email-Enron/email-Enron-nverts.txt",
    "email-Enron/email-Enron-simplices.txt")
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
    #= Returns clique that is weighted =#
    hypergraph = getHypergraph_Weighted()
    clique = hype2Clique_Weighted(hypergraph)
    A1 = cliqueToAdjMatrix_Weighted(clique)
    A = Matrix(A1)
    return A
end


function getHarmonicCentrality(network::Array{Int64, 2}, isWeighted::Bool)
    #= Returns harmonic centrality score for each node as a dictionary =#
    #= Specifiying true to satisfy the isWeighted parameter will calculate
        the harmonic centrality of the network with edge weights being
        considered =#
    G = nx.Graph(network)
    if isWeighted == true
        harmCent=nx.harmonic_centrality(G, distance = "weight")
    else
        harmCent=nx.harmonic_centrality(G)
    end
    return harmCent
end


function getClosenessCentrality(network::Array{Int64,2})
    #= Returns closeness centrality score for each node as a dictionary =#
    G = nx.Graph(network)
    closeCent=nx.closeness_centrality(G)
    return closeCent
end
