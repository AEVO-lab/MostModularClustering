#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#       M O S T    M O D U L A R    C L U S T E R
#
#                       version 2.5
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

import asymmetree.treeevolve   as te
import asymmetree.best_matches as bm
from operator import itemgetter
import numpy    as np
import networkx as nx
import random as rd
import pickle as pkl
import copy
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
from itertools import chain, combinations,product


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#          INPUT-OUTPUT MANAGEMENT FUNCTIONS
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#...........................................................
# This function writes the input files for the HyPPO suite
#...........................................................
def write_inter_cluster_info(scenario_id,cluster,bbh,scenario,matrix):
    #write cluster file
    n1 = scenario_id + ".clusters"
    f1 = open(n1,"w")
    for cc in nx.connected_components(cluster):
        f1.write('<GROUP>\n')
        for node in cc:
            line = str(node) + '_' + str(cluster.nodes[node]['color']) + '\n'
            f1.write(line)
        f1.write('</GROUP>\n')
    f1.close()

    #write score file
    n2 = scenario_id + ".edgelist"
    f2 = open(n2,"w")
    leaves = scenario.genes
    n_matrix = (matrix - np.min(matrix)) / (np.max(matrix) - np.min(matrix)) if np.max(matrix) != np.min(matrix) else matrix/np.max(matrix)

    for v in range(len(leaves)):
        for u in range(len(leaves)):
            line = str(leaves[v].ID)+'_'+str(leaves[v].color) + ' ' +str(leaves[u].ID)+'_'+str(leaves[u].color) +' '+ str(1-n_matrix[v,u]) +'\n'
            f2.write(line)
    f2.close()

    #write species tree
    n3 = scenario_id + '_species_tree.nw'
    f3 = open(n3,"w")
    f3.write(scenario.S.to_newick(distance=False))
    f3.close()

    #write Gene Tree
    n4 = scenario_id +'_true_gene_tree.nw'
    f4 = open(n4,"w")
    f4.write(scenario.TGT.to_newick())
    f4.close()


#...........................................................
# This function prints the adjacency list of a graph including
# the color attributes of the nodes
#...........................................................
def print_ad_list(graph,name):
    print('\n For ',name)
    print('\t NODES \t COLOR \t ADJACENT')
    for v in graph.nodes():
        print('\t',v,'\t',graph.nodes[v]['color'],'\t',list(graph.neighbors(v)))


#..............................................................
# Auxiliary function to extract available colors from a given 
# set of nodes
#..............................................................
def get_species(graph,node_subset):
    colorset = set()
    for n in node_subset:
        if graph.nodes[n]['color'] not in colorset:
            colorset.add(graph.nodes[n]['color'])
    return(list(colorset))



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#       G R A P H    B U I L D E R S
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            
#..............................................................
#  Build best-neighbor graph from scenario
# Takes an asymmetree scenario object and builds the BN graph.
#..............................................................
def bn_from_scenario_matrix(scenario,threshold,matrix,dparam=None,reciprocal=False):
    leaves = scenario.genes
    graph  = nx.DiGraph()
    colors = set()
    
    #normalized matrix for inter-cluster procedure
    n_matrix = (matrix - np.min(matrix)) / (np.max(matrix) - np.min(matrix)) if np.max(matrix) != np.min(matrix) else matrix/np.max(matrix)

    #Add nodes with attributes
    for v in leaves:
        graph.add_node(v.ID,label=v.label,color=v.color)
        colors.add(v.color)

    #Add edges according to similarity threshold (1 - norm_distance)
    for u in range(len(leaves)):
        c_neighbors = {color:[] for color in colors}
        row = matrix[u]
        
        #Add elements to dict
        for v in range(len(leaves)):
            c_neighbors[leaves[v].color].append(v)

        #Normalization by color step
        for c in colors:
            if(c != leaves[u].color):
                c_indices = list(c_neighbors[c])
                m_u = np.take(row,c_indices)
                # If there is only one gene per color, no need to normalize
                if (len(c_indices) == 1):
                    graph.add_edge(leaves[u].ID,leaves[c_indices[0]].ID,similarity= 1.0,score=1-n_matrix[u,c_indices[0]])
                else:
                    n_m = (m_u - min(m_u)) / (max(m_u) - min(m_u)) if max(m_u) != min(m_u) else m_u/max(m_u)

                    for value in range(len(c_indices)):
                        # distance to similarity conversion
                        sim = 1-n_m[value]
                        if( sim >= threshold):
                            graph.add_edge(leaves[u].ID,leaves[c_indices[value]].ID,similarity=sim,score=1-n_matrix[u,c_indices[value]])
                            
    # This part decides whether to keep only bidirectional edges or not
    if (reciprocal == False):
        graph = graph.to_undirected()
    else:
        graph = graph.to_undirected(reciprocal=True)
        
    colors = get_species(graph,graph.nodes())
    
    #Adds colored-neighborhoods attribute to each node    
    for node in graph.nodes():
        for s in colors:
            graph.nodes[node][s] = []
                
        for a in graph.neighbors(node):
            graph.nodes[node][graph.nodes[a]['color']].append(a)    

    if(dparam == None):
        return(graph)
    
    else:
        #In case a d param is given, picks the d most similar neighbors per color
        for node in graph.nodes():
            for s in colors:
                s_neighbors = list(graph.nodes[node][s])
                if(len(s_neighbors) > dparam):
                    edges = sorted([(a,graph.edges[node,a]['similarity']) for a in s_neighbors], key=itemgetter(1), reverse=True)
                    for a,sim in edges[dparam:]:
                        graph.remove_edge(node,a)
                        s_neighbors.remove(a)
                        graph.nodes[a][graph.nodes[node]['color']].remove(node)

                    graph.nodes[node][s] = s_neighbors

        return(graph)


    

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#       COST FUNCTION
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def cost_f(graph,setx,cost_type='default'):
    cost_result = 0

    for vertex in setx:
        non_edges = [y for y in setx if graph.has_edge(vertex,y) == False]
        a_colors  = set(get_species(graph,graph.neighbors(vertex)))
        x_colors  = set(get_species(graph,setx))
        out = a_colors.symmetric_difference(x_colors)

        cost_result = cost_result + (len(non_edges)/ 2) + len(out)

    if(cost_type == 'norm'):
        cost_result = cost_result / len(setx) if len(setx) > 0 else 0 
    
    if(cost_type == 'sqrt'):
        SQRT = math.sqrt(len(setx))
        cost_result = cost_result / SQRT if SQRT > 0 else 0
    
    if(cost_type == 'half'):
        half = len(setx)/2
        cost_result = cost_result / half if half > 0 else 0
        
    return(cost_result)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#            FPT  MAIN  ALGORITHM
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def all_subsets(someiterable):
    s = list(someiterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# Returns a list of tuples
def get_colorful_subsets(graph,vertex):
    n_colors = get_species(graph,graph.neighbors(vertex))
    color_subsets = all_subsets(set(n_colors))
    colorful_subsets = []
    
    for subset in color_subsets:
        some_list = []
        
        for color in subset:
            some_list.append(list(graph.nodes[vertex][color]))
            
        for element in product(*some_list):
            colorful_subsets.append(list(element))
            
    return(colorful_subsets)

def get_Yset(graph,vertex,setx,c_complement):
    some_y = set()
    
    if(len(c_complement) == 0):
        return([[]])
    else:
        closed_n = set(graph.neighbors(vertex))
        closed_n.add(vertex)

        for n in setx:
            for v in graph.neighbors(n):
                if v not in closed_n:
                    if graph.nodes[v]['color'] in c_complement:
                        if v not in some_y:
                            some_y.add(v)
        return(some_y)
        
def get_colorful_Yv(graph,yset,c_complement):
    if (len(c_complement) == 0):
        return([[]]) #empty set
    else:
        reference = dict()
        empty_c   = set()
        
        for c in c_complement:
            reference[c] = list([node for node in yset if graph.nodes[node]['color'] == c])
            if (len(reference[c]) == 0):
                empty_c.add(c)
                
        non_empty_ccom = [ color for color in c_complement if color not in empty_c]
        color_subsets = all_subsets(non_empty_ccom)
        colorful_subsets = []
    
        for subset in color_subsets:
            some_list = []
            for color in subset:
                some_list.append(list(reference[color]))
            
            for element in product(*some_list):
                colorful_subsets.append(list(element))
            
        return(colorful_subsets)

def best_cluster_for_v(graph,vertex,cost_function,cost_t):
    colors = get_species(graph,graph.nodes())
    colors.remove(graph.nodes[vertex]['color'])

    #Get all colorful subsets
    Xv_list = get_colorful_subsets(graph,vertex)
    
    for Xv in Xv_list:
        # current colors in Xv
        c_colors = set(get_species(graph, Xv))
        
        # complement of current colors
        c_complement  = set([color for color in colors if color not in c_colors])

        best_c = [sys.maxsize,[]]

        Y = get_Yset(graph,vertex,Xv,c_complement)
        Yv_list = get_colorful_Yv(graph,Y,c_complement)

        for Yv in Yv_list:
            union = Xv + Yv
            cost  = cost_f(graph,union,cost_t)

            if(cost < best_c[0]):
                best_c[0] = cost
                best_c[1] = union

    return(best_c)

def remove_cluster(graph,cluster):
    #Removes vertices in the cluster from neighbor's colored lists
    for v in cluster:
        for n in graph.neighbors(v):
            graph.nodes[n][graph.nodes[v]['color']].remove(v)
    #Removes cluster vertices from the graph        
    graph.remove_nodes_from(cluster)

#.................................................................
#     m a i n    m m c     f u n c t i o n    
#.................................................................           
def most_modular_cluster(graph,cost_t):
    graph_copy = copy.deepcopy(graph)
    colors = get_species(graph,graph.nodes())
    
    cluster_graph = nx.Graph()
    for v in graph_copy.nodes():
        cluster_graph.add_node(v,label=graph_copy.nodes[v]['label'],
                                       color=graph_copy.nodes[v]['color'])
    flag = True

    singletons = list(nx.isolates(graph_copy))
    graph_copy.remove_nodes_from(singletons)

    while(flag == True):
        best_c = [sys.maxsize,[]]
        
        for v in graph_copy.nodes():
            cost,cluster = best_cluster_for_v(graph_copy,v,cost_f,cost_t)
            if(cost < best_c[0]):
                best_c[0] = cost
                best_c[1] = cluster + [v]

        edges = list(combinations(best_c[1],2))
        cluster_graph.add_edges_from(edges)
        remove_cluster(graph_copy,best_c[1])
        if (graph_copy.order() == 0):
            flag = False
                      
    cluster_graph.add_nodes_from(singletons)

    return(cluster_graph)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#   S P E C T R A L    C L U S T E R I N G
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def remove_tree_like(component):
    flag   = True
    while(flag == True):
        candidates = []
        for v in component.nodes():
            if(component.degree(v) == 1):
                candidates.append(v)

        if(len(candidates) == 0):
            flag = False
        else:
            component.remove_nodes_from(candidates)


def spectral_clustering(graph,a_threshold=0.1):
    graph_copy = copy.deepcopy(graph)
    #Initialize cluster graph
    cluster_graph = nx.Graph()
    for n in graph.nodes():
        cluster_graph.add_node(n,label=graph_copy.nodes[n]['label'],
                                       color=graph_copy.nodes[n]['color'])

    singletons = list(nx.isolates(graph_copy))
    graph_copy.remove_nodes_from(singletons)

    #Main partitioning algorithm
    for cc in nx.connected_components(graph_copy):
        component = graph_copy.subgraph(cc).copy()
        remove_tree_like(component)
        component.remove_nodes_from(list(nx.isolates(component)))

        if component.order() > 0:
            c_list = []

            if(nx.is_connected(component) == True):
                c_list.append(component)
            else:
                for subcc in nx.connected_components(component):
                    c_list.append(component.subgraph(subcc).copy())

            flag = True

            while(flag == True):
                sub_component = c_list.pop()
                a_connectivity = nx.algebraic_connectivity(sub_component)/sub_component.order()
            
                if(a_connectivity > a_threshold):
                    cluster_graph.add_edges_from(list(sub_component.edges()))
                else:
                    fiedler_vector = nx.fiedler_vector(sub_component)
                    nodes = np.array(sub_component.nodes())
                    #weak nodal domain
                    c1 = sub_component.subgraph(nodes[np.where(fiedler_vector >= 0)]).copy()
                    c1.remove_nodes_from(list(nx.isolates(c1)))
                    
                    c2 = sub_component.subgraph(nodes[np.where(fiedler_vector < 0)]).copy()
                    c2.remove_nodes_from(list(nx.isolates(c2)))

                    if(c1.order() > 1):
                        if nx.number_connected_components(c1) > 1:
                            for component in nx.connected_components(c1):
                                if c1.subgraph(component).order() > 1:
                                    c_list.append(c1.subgraph(component))
                        else:
                            c_list.append(c1)
                        
                    if(c2.order() > 1):
                        if nx.number_connected_components(c2) > 1:
                            for component in nx.connected_components(c2):
                                if c2.subgraph(component).order() > 1:
                                    c_list.append(c2.subgraph(component))
                        else:
                            c_list.append(c2)

                if(len(c_list) == 0):
                    flag = False

    return(cluster_graph)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#        HyPPO's  OCR program output interpretation
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def graph_editing_w_HyPPO(outputfile,cluster_graph):
    interc_edges_ortho = set()
    interc_edges_para  = set()

    with open(outputfile) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            if (cnt == 2):
                line = line.rstrip("\n")
                parsed = line.replace('RELATIONS=','')
                as_list = parsed.split(';;')

                for x in as_list:
                    p = x.split('\t')
                    gene_1 = p[0]
                    gene_2 = p[1]
                    rel = p[2]
                    p1 = gene_1.split('_')
                    p2 = gene_2.split('_')
                    edge = (int(p1[0]),int(p2[0]))
                    
                    if(rel=='Orthologs'):
                        interc_edges_ortho.add(edge)
                    else:
                        interc_edges_para.add(edge)
            line = fp.readline()
            cnt += 1

    for edge in cluster_graph.edges():
        if edge in interc_edges_para:
            print("\t Paralogous edge ",edge," detected in cluster graph")
            sys.exit()
            cluster_graph.remove_edge(*edge)

    for edge in interc_edges_ortho:
        if (cluster_graph.has_edge(*edge) == False):
            cluster_graph.add_edge(*edge)


            
                
                
                
    
