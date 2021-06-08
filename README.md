# MostModularClustering

The source code and a clean readme will be available by the Recomb-CG deadline.

`Most Modular Clustering` is a library that contains a python implementation of the FPT algorithm **Most Modular Clustering (MMC)** which can also be combined with the [HyPPO](https://github.com/AEVO-lab/HyPPO) suite (MMC+OCR) to predict orthologs in simulated data.

![for_git](https://user-images.githubusercontent.com/20115213/121257682-987c2880-c873-11eb-9366-db013024745f.png)



## Dependencies
`Most Modular Clustering` has the following dependencies:

| package name | version |
|:------------:|:-------:|
|  asymmetree  |  1.0.2  |
|   networkx   |   2.5   |
|    numpy     |  1.20.1 |
|    pandas    |  1.2.3  |
|     scipy    |  1.6.1  |
|    seaborn   |  0.11.1 |
|  matplotlib  |  3.3.4  |

## Usage and Description
To use `Most Modular Clustering` just download `most_modular_cluster_2_5.py` to your current working directory and just import the library as follows: 

```python
import most_modular_cluster_2_5 as mmc
```
To generate a phylogenetic scenario, we used the [`AsymmeTree`](https://github.com/david-schaller/AsymmeTree) library. From this library, it is possible to simulate noisy evolutionary distance data to test our clustering algorithms. 

### (k,d)-Best Neighbor Graph
From `AsymmeTree`'s distance matrix and `scenario` object, it is possible to construct a (k,d)-Best Neighbor Graph where user should specify `threshold` (between 0.0 and 1.0) and optionally:
* `dparam` : max. number of neighbors from each color, otherwise all edges above threshold are taken. 
* `reciprocal`:`True` saves only the reciprocal edges, i.e. if gene a is above similarity threshold with b and viceversa, ab will be an edge in the final (k,d) graph.

```python
bn = mmc.bn_from_scenario_matrix(scenario,threshold,matrix,dparam=d,reciprocal=False)
```
This graph will serve as input for our clustering algorithms.

### Most Modular Clustering Algorithm
`Most Modular Clustering`assumes input graph is a `networkx` undirected graph where every node has a `color` attribute
and every node has a list of neighbors corresponding to each available color. This function fill output a `networkx` undirected graph where every connected component is a clique.
 
```python
mmc_graph = mmc.most_modular_cluster(graph,cost_function)
```
The user should specify the cost function to be used by the algorithm. Some of the available cost functions are:
|  function | description                                                   |
|:---------:|---------------------------------------------------------------|
| `default` | default cost as described in (manuscript)                     |
|   `norm`  | default cost divided by the size of candidate set             |
|   `sqrt`  | default cost divided by the square root of candidate set size |
|   `half`  | default cost divided by half of the candidate set size        |

### MMC + OCR
It is possible to combine our algorithm with the **Orthology Cluster Recovery (OCR)** module of the `HyPPO` suite which computes the missing orthology relations that could exist between clusters. To output the necessary files for OCR use:

```python
write_inter_cluster_info(scenario_id,mmc_graph,bn,scenario,matrix)
```
where `scenario_id` is an identifier for the files and `mmc_graph` should be the cluster graph obtained in the previous step. When the files are created the user can run OCR on them. Afterwards, it is possible to add the missing edges directly to the `mmc_graph` by using:

```python
graph_editing_w_HyPPO(outputfile,mmc_graph)
```
where `outputfile` is the file name given to OCR's output.

### Spectral Clustering Algorithm
This algorithm is based on the clustering part of [ProteinOrtho](https://www.bioinf.uni-leipzig.de/Software/proteinortho/), it decomposes  the  input  graph  until  every  connected  component  reachesan  algebraic  connectivity  above  a  certain  threshold `a_threshold` which is an optional parameter. It returns an undirected graph with the resulting clusters.
```python
spc_graph = spectral_clustering(graph,a_threshold=0.1)
```
