# HyperharmonicCentrality
 The repository for our ESMI 2020 Research on calculating Harmonic Centrality on Higher Order Sets

# Summary <h1>
In this project contains methods for converting data files containing
hypergraphs into cliques and then transforming the cliques into their
adjacency matrix representation. We then model several different networks with
their centralization visualized. Finally, we run a couple of statistic calculations
to analyze the network's different type of centrality measures.

## Visuals <h2>
**Graphical Results on the Enron Email Dataset are as follows:**

*Unweighted Harmonic Centrality Visualization*
![(UnWeighted) Harmonic Centrality](/GraphVisuals/CliqueHarmonicRED.png)

*Weighted Harmonic Centrality Visualization*
![(Weighted) Harmonic Centrality](/GraphVisuals/CliqueHarmonic_WeightedRED.png)

*Unweighted Closeness Centrality Visualization*
![(UnWeighted) Closeness Centrality](/GraphVisuals/CliqueClosnessRED.png)

## Statistics <h2>
**Analysis Summary on the Enron Email Dataset:**

* *Nodes, Edges, and Hyperedges*
    * Nodes      - 143
    * Edges      - 1800
    * Hyperedges - 10551

* *Euclidean Normalization (0-2 scale)*
    * Unweighted and Weighted Harmonic Centrality                        - 0.99899
    * Unweighted Harmonic Centrality and Unweighted Closeness Centrality - 0.1869

*Spearman Correlation Coefficient on Unweighted and Weighted Harmonic Rankings*
![SpearmanCorrCoeff](/GraphVisuals/6241-spearman-HarmonicRankings-rankcorr.png)
