using DelimitedFiles
using Random
using  PyCall
nx=pyimport("networkx")

function hyp2Cliq(dataFile1::String, dataFile2::String)
    #=
    Returns: clique representing the projected graph of
                a hypergraph

    Parameter: File containing hypergraph
<<<<<<< HEAD
    Precondition: files exist and are in the same folder, first file lists the
    hyperedges and second file lists the nodes
=======
    Precondition: Files exist and ar in the same folder
>>>>>>> 20e388b0a2c945b422a382d0ccede870c15e40df
    =#
    #FileReading Hypergraph (if we have two Files)

    #first file's line number* is its position in the matrix
    #the number on that line* is how many nodes are in the edge
    #go into second file and take that many nodes from last (keep track of
    #position last taken from)

    #final array like this:
    # hyper = [ [nodes in hyperedge1]; [nodes in hyperedge2]; ... ;[nodes in hyperedge n] ]

    edges = readdlm(dataFile1)
<<<<<<< HEAD
    nodes = readdlm(dataFile2)
=======
    nodes = readdlm(datafile2)
>>>>>>> 20e388b0a2c945b422a382d0ccede870c15e40df
    numedges = size(edges,1) #number of hyperedges / rows in hypergraph matrix
    numnodes = size(nodes,1) #total number of hypernodes


<<<<<<< HEAD
    #initialize hypergraph as an array of arrays
    hyper = hyper=[[1,2]]

    #go through each line
    i = 1 #line pos we're at in second file
    for line in 1:numedges
        take = edges[line] #number of nodes to take from other file
        Edge = nodes[i:(i-1)+Int64(take)] #makes an array representing hyperedge with its nodes
        i += Int64(take) #update line postion in second file
        push!(hyper, Edge) #place hyperedge array in overall hypergraph
    end
    hyper=hyper[2:end]
=======
    #first crate a hypergraph array and zero it out
    hyper = zeros(Int64, numedges ,1)

    #go through each line

    i = 1 #line pos we're at in second file
    for line in 1:numedges
        take = edges[line] #number of nodes to take from other file
        Edge = nodes[i:i+take] #makes an array representing hyperedge with its nodes
        i += take #update line postion in second file
        hyper[line] = Edge #place hyperedge array in overall hypergraph
>>>>>>> 20e388b0a2c945b422a382d0ccede870c15e40df


    # Hypergraph (Converting if hyperedge is line by line)
    # line2vec(line) = parse.(Int, split(line, '\t' ))
    # hyper=[line2vec(line) for line in eachline(dataFile;
    #         keep = false) if !isempty(line)]


#     # Conversion to clique
#     clique = Vector{Array{Int64, 1}}
#
#
#
#     n = size(hyper, 1)
#     for k = 1:n
#         m = size(hyper, 2)
#         temp=zeros(m)
#         for l=1:m
#             temp[l]
#     end
#
#
#     # Once we have clique, this makes the adjency matrix of clique
#     A = max.(clique, clique')
#
 end


hyp2Cliq("email-Enron/email-Enron-nverts.txt", "email-Enron/email-Enron-simplices.txt")
