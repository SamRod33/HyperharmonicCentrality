using DelimitedFiles
using Random
using  PyCall
nx=pyimport("networkx")

function hyp2Cliq(dataFile::String)
    #=
    Returns: clique representing the projected graph of
                a hypergraph

    Parameter: File containing hypergraph
    Precondition: File that contains a set of nodes
                    in the Hyperedge on one line,
                    and each line representing a hyperedge
    =#

    # Hypergraph
    line2vec(line) = parse.(Int, split(line, '\t' ))
    hyper=[line2vec(line) for line in eachline(dataFile;
            keep = false) if !isempty(line)]
    # Conversion to clique
    clique=Int64[]

    # Once we have clique, this makes the adjency matrix of clique
    n = # Fill in
    data = sparse(convert(Vector{Int64}, data[:,1]),
               convert(Vector{Int64}, data[:,2]), 1, n, n)
    data = max.(data, data')

end
