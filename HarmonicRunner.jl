using DelimitedFiles
using Random
using  PyCall
nx=pyimport("networkx")

function hyp2Cliq(dataFile::String, dataFile::String)
    #=
    Returns: clique representing the projected graph of
                a hypergraph

    Parameter: File containing hypergraph
    Precondition:
    =#

    #FileReading Hypergraph (if we have two Files)

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
