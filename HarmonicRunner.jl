using DelimitedFiles
using Random
using  PyCall
nx=pyimport("networkx")

function hyp2Cliq(dataFile1::String, dataFile2::String)
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
        Edge = nodes[i:(i-1)+Int64(take)] #makes an array representing hyperedge with its nodes
        i += Int64(take) #update line postion in second file
        push!(hyper, Edge) #place hyperedge array in overall hypergraph
    end
    hyper=hyper[2:end]


    # Hypergraph (Converting if hyperedge is line by line)
    # line2vec(line) = parse.(Int, split(line, '\t' ))
    # hyper=[line2vec(line) for line in eachline(dataFile;
    #         keep = false) if !isempty(line)]


    # Conversion to clique
    clique = Array[]

    n = size(hyper, 1)
    for hyperedge = 1:n
        m = size(hyper[hyperedge,:], 2)
        for currentNode=1:m
            for otherNodes=currentNode+1:m
                i = hyper[hyperedge][1,currentNode]
                j= hyper[hyperedge][1,otherNode]
                push!(clique, [i j])
            end
        end
    end

    return clique
    # Once we have clique, this makes the adjency matrix of clique
 end


hyp2Cliq("email-Enron/email-Enron-nverts.txt", "email-Enron/email-Enron-simplices.txt")
