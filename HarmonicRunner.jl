using DelimitedFiles
using Random
using  PyCall
nx=pyimport("networkx")

function hyp2Cliq(dataFile::String, dataFile::String)
    #=
    Returns: clique representing the projected graph of
                a hypergraph

    Parameter: File containing hypergraph
    Precondition: Two
    =#

    #FileReading Hypergraph (if we have to Files)

    # Hypergraph (Converting if hyperedge is line by line)
    # line2vec(line) = parse.(Int, split(line, '\t' ))
    # hyper=[line2vec(line) for line in eachline(dataFile;
    #         keep = false) if !isempty(line)]

    # Conversion to clique
    clique=Int64[]

    # Once we have clique, this makes the adjency matrix of clique
    A = max.(clique, clique')

end
