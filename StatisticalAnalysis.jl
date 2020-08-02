#=
Module for analyzing the Enron Email Dataset

This module runs several statistical calculations, including the L_2 Norm (or
Euclidean Norm) between:
                1) Unweighted Harmonic Centrality and Closeness Centrality
                2) Weighted Harmonic Centrality and Weighted Harmonic Centrality

Additionally, we analyze the Unweighted and Weighted Harmonic Centrality node
rankings by calculating its Spearman Correlation Coefficient

Author: Samuel Rodriguez (sar325@cornell.edu)
Date:   08.01.20
=#

using PyPlot: savefig, clf, semilogx, gca, legend, xlabel, ylabel, title
include("HyperNetwork.jl")

#######################################################################
#------------------DATA for Enron Data Set:---------------------------#
#######################################################################
# nodes, edges, hyperedges
# nodes      - 143
# edges      - 1800
# hyperedges - 10551




#######################################################################
#---------------------------------------------------------------------#
#---------------------STATISTICS ANALYSIS-----------------------------#
#---------------------------------------------------------------------#
#######################################################################

#######################################################################
#-----Calculating L_2 Norm for harmonic and closeness centrality------#
#-----(How close are the two datasets)--------------------------------#
#######################################################################
hCent = getHarmonicCentrality(getAdjMatrix(), false)
ordered_hCent=sort(hCent)
hCentScores = [t[2] for t in ordered_hCent]
hCentScoresNorm = hCentScores / maximum(hCentScores)

cCent = getClosenessCentrality(getAdjMatrix())
ordered_cCent=sort(cCent)
cCentScores = [t[2] for t in ordered_cCent]
cCentScoresNorm = cCentScores/ maximum(cCentScores)

norm(hCentScoresNorm-cCentScoresNorm,2)
### L_2 Norm --> 0.1869 result on a normalized scale (i.e from 0-2 scale) ###


##############################################################################
######-------------- WHAT IF WE DIDNT MAKE ALL WEIGHTS ONE ------------#######
##############################################################################
hCentW=getHarmonicCentrality(getAdjMatrix_Weighted(), true)
ordered_hCentW=sort(hCentW)
hCentWScores = [t[2] for t in ordered_hCentW]
hCentWScoresNorm = hCentWScores / maximum(hCentWScores)

hCentNotW=getHarmonicCentrality(getAdjMatrix(), false)
ordered_hCentNotW=sort(hCentNotW)
hCentNotWScores = [t[2] for t in ordered_hCentNotW]
hCentNotWScoresNorm = hCentNotWScores / maximum(hCentNotWScores)

norm(hCentNotWScoresNorm-hCentWScoresNorm,2)
# L_2 Norm --> 0.99899 (0 - 2 range, greater distance from 0
#                       indicates two data sets are more different)



########################################################################
#---- Spearman Correlation Coefficient (compares ranking of nodes) ----#
#---- between unweighted and weighted harmonic centralities -----------#
########################################################################

function rank_corrbetweenness(name::Any,dataset::Any, index_start::Int64)
    # Compares ranking of nodes between two sets of data
    cec_c=dataset[:,1]
    n=length(cec_c)
    cec_c[cec_c]=[i for i in 1:n]
    hec_c=dataset[:,2]
    hec_c[hec_c]=[i for i in 1:n]
    start = log10(index_start - 1)
    finish = log10(length(cec_c) - 1)
    ran = [convert(Int64, round(v)) for v in 10 .^ range(start, stop=finish, length=500)]
    sp = sortperm(cec_c)
    cec_c, hec_c= cec_c[sp], hec_c[sp]
    cec_hec = [corspearman(cec_c[(end-s):end], hec_c[(end-s):end]) for s in ran]
    sp = sortperm(hec_c)
    cec_c, hec_c = cec_c[sp], hec_c[sp]
    hec_cec = [corspearman(hec_c[(end-s):end], cec_c[(end-s):end]) for s in ran]
    clf()
    xs = collect(ran) .+ 1
    semilogx(xs, cec_hec, linestyle="-",  lw=1,    label="with respect to weighted")
    semilogx(xs, hec_cec, linestyle=":",  lw=2.25, label="with respect to unweighted")
    fsz = 20
    ax = gca()
    ax.tick_params(axis="both", labelsize=fsz-2, length=7, width=1.5)
    ax.tick_params(axis="x",    which="minor", length=4, width=1)
    legend(fontsize=fsz-10, frameon=false)
    title("$name", fontsize=fsz)
    xlabel("Number of top-ranked elements", fontsize=fsz)
    ylabel("Spearman’s rank corr. coeff.", fontsize=fsz)
    savefig("6241-spearman-$name-rankcorr.png", bbox_inches="tight")
end


# Creates lists by ordering nodes from highest centrality score to lowest
hCentWRank = [n[1] for n in sort(hCentW, byvalue=true, rev=true)]
hCentNotWRank = [n[1] for n in sort(hCentNotW, byvalue=true, rev=true)]

# Combines both lists where each list is put in its own column
ranks = [ [(n+1) for n in hCentWRank] [(k+1) for k in hCentNotWRank] ]

# Calculates Spearman Correlation Coefficient and Plots it
rank_corrbetweenness("HarmonicRankings", ranks, 10)
