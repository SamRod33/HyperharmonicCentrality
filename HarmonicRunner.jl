using DelimitedFiles
using Random
using  PyCall
nx=pyimport("networkx")

function hyp2Cliq(dataFile1::String, dataFile2::String)
    #=
    Returns: clique representing the projected graph of
                a hypergraph

    Parameter: File containing hypergraph
    Precondition:
    =#
    #FileReading Hypergraph (if we have two Files)
    
    #first file's line number* is its position in the matrix
    #the number on that line* is how many nodes are in the edge
    #go into second file and take that many nodes from last (keep track of
    #position last taken from)

    #final matrix like this:
    # hyper = [ [nodes in hyperedge1]; [nodes in hyperedge2]; ... ; [nodes in hyperedge n] ]

    edges = readdlm(dataFile1)
    nodes = readdlm(datafile2)
    rows = size(nodes,1) #number of hyperedges / rows in hypergraph matrix

    for line in 1:rows
        nodes = edges[line] #number of nodes to take from other file



    #first crate a matrix and zero it out?
    hyper = zeros([T=Float64,] Int64::2,3)

    hyper # Hypergraph Variable

    # Hypergraph (Converting if hyperedge is line by line)
    # line2vec(line) = parse.(Int, split(line, '\t' ))
    # hyper=[line2vec(line) for line in eachline(dataFile;
    #         keep = false) if !isempty(line)]

    # Conversion to clique
    clique = Vector{Array{Int64, 1}}



    n = size(hyper, 1)
    for k = 1:n
        m = size(hyper, 2)
        temp=zeros(m)
        for l=1:m
            temp[l]
    end


    # Once we have clique, this makes the adjency matrix of clique
    A = max.(clique, clique')

end
