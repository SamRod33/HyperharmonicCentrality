#=
Module for displaying network

This module provides a couple of methods for visualing a network given
its adjecency matrix representation. The network's Harmonic, Closeness, and
Weighted Harmonic Centralities visualizations are supported.

Author: Samuel Rodriguez (sar325@cornell.edu)
Date:   08.01.20
=#

import Pkg; Pkg.add("GraphPlot")
using GraphPlot
using LightGraphs
using Colors
using Cairo
using Compose: draw, cm, PNG
include("HyperNetwork.jl")


#############################################################################
#-------------------- GRAPHING NETWORKS ------------------------------------#
#############################################################################
function graphHarmonicNetwork_Weighted(network::Array{Int64, 2})
    #=
    Returns: graph of clique with node color scaled by
            their *weighted* harmonic centrality score
    =#
    Random.seed!(12345678)  # plot same each time
    harmCent = getHarmonicCentrality(network, true)
    n = length(keys(harmCent))-1
    cent=[harmCent[i] for i in 0:n]
    cent = cent / maximum(cent)
    g = SimpleGraph(network)
    nodefill = [RGBA(i,0.0,0.0,i) for i in cent]
    layout=(args...)->spring_layout(args...; C=35)
    p= gplot(g, layout=layout, nodefillc=nodefill,
            edgestrokec=colorant"grey", nodelabel = 1:143,
            NODELABELSIZE=3, NODESIZE=0.040)
    draw(PNG("CliqueHarmonic_WeightedRED.png", 16cm, 16cm), p)
end

function graphHarmonicNetwork(network::Array{Int64, 2})
    #=
    Returns: graph of clique with node color scaled by
            their harmonic centrality score
    =#
    Random.seed!(12345678)  # plot same each time
    harmCent = getHarmonicCentrality(network, false)
    n = length(keys(harmCent))-1
    cent=[harmCent[i] for i in 0:n]
    cent = cent / maximum(cent)
    g = SimpleGraph(network)
    nodefill = [RGBA(i,0.0,0.0,i) for i in cent]
    layout=(args...)->spring_layout(args...; C=35)
    p= gplot(g, layout=layout, nodefillc=nodefill,
            edgestrokec=colorant"grey", nodelabel = 1:143,
            NODELABELSIZE=3, NODESIZE=0.040)
    draw(PNG("CliqueHarmonicRED.png", 16cm, 16cm), p)
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
graphHarmonicNetwork_Weighted(getAdjMatrix_Weighted())
