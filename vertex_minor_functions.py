#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scipy
import networkx as nx #A python library that is helpful when dealing with graphs 
import matplotlib.pyplot as plt


# In[2]:


## some useful functions 

##draw the grid graph, can choose to lock the vertices into grid structure 
def draw_graph( graph, position = False, plot_size = (5,5)):

    plt.figure(figsize= plot_size)

    if position == False: 

        nx.draw(
        graph,
        with_labels=True,
        node_color="skyblue",
        node_size=600,
        font_size=8,
        font_weight = 'bold',
        edge_color="gray",)

    else: 
        nx.draw(
            graph,
            pos= position,
            with_labels=True,
            node_color="skyblue",
            node_size=600,
            font_size=8,
            font_weight = 'bold',
            edge_color="gray",)      

    plt.show()

#function to determine where each node needs to be for gridlock for a grid defined by integers
def convert_position_to_integers (integer_Graph, cols ):
    pos_int = {}
    for node in integer_Graph.nodes():
        node_mod = node- 1
        col_idx = node_mod % cols
        row_idx = node_mod // cols
        pos_int[node] = (col_idx, -row_idx)
    return pos_int


# In[3]:


##perform local complementation on a node using a graph defined by integers
def local_complementation_integers (vertex, Graph, show_graph = False, gridlock = False, plot_size = (3,3), position_int = None):    
    # copy the graph 
    G = Graph.copy()

    #extract neighbor subgraph of vertex 
    neighbors = list(G.neighbors(vertex))
    G_sub = G.subgraph(neighbors)

    #complement the subgraph 
    G_sub_complement = nx.complement(G_sub)

    #delete existing edges between neighbors 
    G.remove_edges_from(G_sub.edges())

    #add complement edges between naeighbors
    G.add_edges_from(G_sub_complement.edges())

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G    


# In[4]:


#Functions for performing Graph Operations/Pauli Measurements
##Each function has the option to draw the resulting graph after performing the operation

def Z_graph_operation(vertex_to_delete, graph_int, show_graph = False, gridlock = False, 
                      plot_size = (3,3), position_int = None): 

    G = graph_int.copy()

    #delete vertex
    G.remove_node(vertex_to_delete)

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G


def Y_graph_operation(vertex_to_delete, graph_int, show_graph = False, gridlock = False, 
                      plot_size = (3,3), position_int = None): 

    G = graph_int.copy()

    #perform LC on vertex 
    G = local_complementation_integers(vertex_to_delete, G)

    #delete vertex
    G.remove_node(vertex_to_delete)

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G


def X_graph_operation(vertex_to_delete, graph_int, pivot = False, show_graph = False, gridlock = False, 
                      plot_size = (3,3), position_int = None): 

    neighbors = list(graph_int.neighbors(vertex_to_delete))

    ##check to see if the vertex has any neighbors 
    if not neighbors:
        print ("vertex has no neighbors, unable to perform operation")
        return

    #pick a pivot if unspecified 
    if pivot == False:
        pivot = neighbors[0]

    #ensure the pivot choice is a neighbor to the desired vertex 
    if pivot not in neighbors:
        print("pivot is not a neighbor of vertex, unable to perform operation")
        return graph_int

    G = graph_int.copy()

    #perform LC's on vertex and pivot  
    G = local_complementation_integers(pivot, G)
    G = local_complementation_integers(vertex_to_delete, G)
    G = local_complementation_integers(pivot, G)

    #delete vertex
    G.remove_node(vertex_to_delete)

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G


# In[16]:


#functions for protocols on 3xn and 2xn grids 

#for a 3xn graph 
def crossing_protocol(vertex_a, graph_int, cols, pivot = False, show_graph = False, gridlock = False, 
                      plot_size = (3,3), position_int = None): 

    neighbors = list(graph_int.neighbors(vertex_a))

    #check to see if vertex a is a valid choice
    required_neighbors = [vertex_a - cols, vertex_a -1, vertex_a + 1, vertex_a + cols]
    for neighbor in required_neighbors:
        if neighbor not in neighbors:
            print ("vertex a is not not a valid option, please pick another vertex.")
            return

    G = graph_int.copy()

    #perform X measurement on vertex a using a-n as the pivot 
    G= X_graph_operation(vertex_a, graph_int, pivot = vertex_a-cols, show_graph= True, gridlock = True, plot_size= (10,3), position_int= pos_int)

    #perform vertex deletions 
    G = Z_graph_operation( vertex_a-cols,G )
    G = Z_graph_operation( vertex_a+cols,G )

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G


# for a 3xn graph
def flipping_protocol(graph_int, cols, show_graph = False, gridlock = False, 
                      plot_size = (3,3), position_int = None): 

    #check to see if grid is sufficiently large 
    if cols < 4:
        print ("Grid is not large enough.")
        return

    G = graph_int.copy()

    #perform X and Z measurements 
    G =  X_graph_operation( cols -1, G, pivot= cols )
    G =  X_graph_operation( 2*cols -1, G, pivot= 3*cols-1 )
    G = Z_graph_operation( 3*cols-1 ,G )

    #perform LC's
    G = local_complementation_integers( 2*cols ,G )
    G = local_complementation_integers( 3*cols ,G )
    G = local_complementation_integers( 2*cols ,G )

    #perform Y measurements 
    G = Y_graph_operation(cols-2 ,G )
    G = Y_graph_operation(2*cols-2 ,G )
    G = Y_graph_operation(3*cols-2 ,G )

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G

# for a 3xn graph
def extract_corners(graph_int, cols, show_graph = False, gridlock = False, 
                      plot_size = (3,3), position_int = None, invert = False): 

    #check to see if grid is sufficiently large 
    if cols < 3:
        print ("Grid is not large enough.")
        return

    G = graph_int.copy()

    #perform X Y Z measurements for all n
    if invert == True: 
        G =  Y_graph_operation( 2*cols+2, G)
        G =  X_graph_operation( cols +2, G, pivot= 2 )
        G = Z_graph_operation( 2 ,G )
        G = Y_graph_operation( cols+3 ,G )
        G = Y_graph_operation( 3 ,G )
    else:
        G =  Y_graph_operation( 2, G)
        G =  X_graph_operation( cols +2, G, pivot= 2*cols+2 )
        G = Z_graph_operation( 2*cols+2 ,G )
        G = Y_graph_operation( cols+3 ,G )
        G = Y_graph_operation( 2*cols+3 ,G )

    #if n=3 the protocol is over 
    if cols == 3: 
            #plot graph 
        if show_graph== True:
            if gridlock == False:
                draw_graph(G, plot_size = plot_size)
            else:
                draw_graph(G, position_int, plot_size = plot_size)

        return G


    #if n>=4 finish the protocol 
    if invert == True:
        G = Z_graph_operation(cols+4 ,G )
        G = Z_graph_operation(4 ,G )
    else:
        G = Z_graph_operation(cols+4 ,G )
        G = Z_graph_operation(2*cols+4 ,G )

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G


# In[5]:


#perform the X Protocol on a grid of integers to extract a bell pair. Must specify the desired ordered repeater line
def x_protocol_grid_integers(graph_int, ordered_repeater_line, bell_pair, show_graph = False, gridlock = False, plot_size = (3,3),
                             position_int = None, show_X_measurement = False, gridlock_X = False, show_Z_measurement = False,
                             gridlock_Z = False):

    G = graph_int.copy()

    #itratively perform x measurements along the repeater line 
    for i in range (len(ordered_repeater_line)-2):
        #pivots are specified using the following: let V_i be the ith vertex in the repeater line. Then if X(V_1), then pivot = V_0; 
        # if x(v_{n-1}), then pivot = V_n; for any intermeidary, use an adjacent pivot on the repeater line. 

        if i == 0: 
            p = ordered_repeater_line[0]
        elif i == len(ordered_repeater_line)-3:
            p = ordered_repeater_line[-1]
        else :
            neighbors = list(G.neighbors(ordered_repeater_line[i+1]))
            pivots = [v for v in neighbors if v in ordered_repeater_line]

            if not pivots:
                raise ValueError(f"No repeater neighbor found for vertex {ordered_repeater_line[i+1]}")

            p = pivots[0]


    G = X_graph_operation(ordered_repeater_line[i+1], G, pivot = p, show_graph= show_X_measurement, gridlock= gridlock_X, 
                          plot_size= plot_size)

    #find all the neighbors of the bell pair 
    neighboring_vertices = list( G.neighbors(bell_pair[0])) + list( G.neighbors(bell_pair[1])) 
    neighboring_vertices = [v for v in neighboring_vertices if v not in bell_pair]

    #perform z measurements on all neighbors
    for v in neighboring_vertices: 
        G = Z_graph_operation(v, G, show_graph = show_Z_measurement, gridlock= gridlock_Z, plot_size= plot_size)

    #plot graph 
    if show_graph== True:
        if gridlock == False:
            draw_graph(G, plot_size = plot_size)
        else:
            draw_graph(G, position_int, plot_size = plot_size)

    return G

