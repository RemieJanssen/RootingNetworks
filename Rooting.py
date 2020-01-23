import networkx as nx
import itertools
import matplotlib.pyplot as plt
import sys
from copy import deepcopy



###############################################
###############    CLASSES   ##################
###############################################
# Functions that check whether a network is in a specific rooted class

#Checks if a directed graph is a rooted tree. If so, it returns the root node
def IsRootedTree(digraph):
    if len(digraph.edges)!=len(digraph.nodes)-1 or not nx.is_weakly_connected(digraph):
        return False
    rootFound = False
    root = False
    for node in digraph:
        if digraph.in_degree(node)==0:
            if rootFound:
                return False
            rootFound = True
            root = node
        if digraph.in_degree(node)>1:
            return False
    return root


#A trivial function, that returns True always
#In particular, for each phylogenetic network, it returns True.
def ClassAllNetworks(network):
    return True
    
#A function that returns True iff the given network is tree-child.    
def ClassTreeChild(network):
    for node in network.nodes:
        tc_node = False
        if network.out_degree(node) == 0:
            tc_node = True
        for child in network.successors(node):
            outdegChild = network.out_degree(child)
            if outdegChild == 0 or outdegChild == 2:
                tc_node = True
        if not tc_node:
            return False
    return True

#A function that returns True iff the given network is stack-free.    
def ClassStackFree(network):
    for node in network.nodes:
        sf_node = True
        if network.out_degree(node) == 1:
            for child in network.successors(node):
                outdegChild = network.out_degree(child)
                if outdegChild == 1:
                    sf_node = False
        if not sf_node:
            return False
    return True
    
#A function that returns True iff the given network is tree-based.    
def ClassTreeBased(network):
    for node in network.nodes:
        if network.out_degree(node)==1:
            if AtWFence(network,node):
                return False
    return True
    

#A function that returns True iff the given network is orchard.    
def ClassOrchard(network):
    redNw = deepcopy(network)
    leaves = set()
    for node in redNw:
        if redNw.out_degree(node)==0:
            leaves.add(node)
        if redNw.in_degree(node)==0:
            root = node        
    redNw.add_edge(-2,root)
    done = False
    while not done:
        checkedAllLeavesMaybe = True
        for l in leaves:
            pair = IsSecondInPair(redNw,l)
            if pair:
                reduced = ReducePair(redNw,*pair)
                if reduced == "C":
                    leaves.remove(pair[0])
                checkedAllLeavesMaybe = False
                break
        if len(redNw.edges)==1:
            return True
        done = checkedAllLeavesMaybe
    return False    


###############  Some subroutines for the class checkers   ###################

#Checks if the network has a W-fence with the given node at one of the endpoints of the fence
def AtWFence(network,reticulation):
    for child in network.successors(reticulation):
        currentNode = child
    previousNode = reticulation
    currentlyAtTop=False
    done = False
    while not done:
        currentOutDegree = network.out_degree(currentNode)
        if currentOutDegree==0:
            return False
        if currentOutDegree==1:
            if currentlyAtTop:
                return True
            for node in network.predecessors(currentNode):
                if node!=previousNode:
                    nextNode = node
        if currentOutDegree==2:
            if not currentlyAtTop:
                return False
            for node in network.successors(currentNode):
                if node!=previousNode:
                    nextNode = node
        previousNode, currentNode = currentNode, nextNode
        currentlyAtTop = not currentlyAtTop
    return False, "Error"


#Reduces the pair (x,y) in network if it is reducible in network.
def ReducePair(network,x,y):
    for node in network.predecessors(x):
        px = node
    for node in network.predecessors(y):
        py = node
    if px == py:
        for node in network.predecessors(py):
            ppy = node
        network.remove_edges_from([(ppy,py),(py,x),(py,y)])
        network.add_edge(ppy,y)
        return "C"
    if py in network.predecessors(px):
        for node in network.predecessors(py):
            ppy = node
        for node in network.predecessors(px):
            if node != py:
                ppx = node
        network.remove_edges_from([(ppy,py),(py,px),(py,y),(ppx,px),(px,x)])
        network.add_edges_from([(ppy,y),(ppx,x)])
        return "RC"
    return False


# Returns a reducible pair of network with x as second element if it exists, returns False otherwise.
def IsSecondInPair(network,x):
    for node in network.predecessors(x):
        px = node
    cpx = False
    for node in network.successors(px):
        if node != x:
            cpx = node
    if type(cpx)==bool and not cpx:
        return False
    if network.out_degree(cpx)==0:
        return (cpx,x)
    if network.out_degree(cpx)==1:
        for node in network.successors(cpx):
            if network.out_degree(node)==0:
                return (node,x)
    return False



###############################################
############  BASIC ORIENTATION  ##############
###############################################



def OrientationAlgorithmBinary(network,rootEdge,reticulations):
    numberNetworkEdges = len(network.edges)
    network.remove_edges_from([rootEdge])
    if len(network.nodes)+len(reticulations)!=numberNetworkEdges+1:
        network.add_edges_from([rootEdge])
        return False
    diNetwork=nx.DiGraph()
    diNetwork.add_nodes_from(network.nodes)
    diNetwork.add_node(-1)
    diNetwork.add_edges_from([(-1,rootEdge[0]),(-1,rootEdge[1])])
    readyNodes = (network.nodes - reticulations) & set(rootEdge)
    numberOrientedEdges = 2
    while len(readyNodes)>0: #subdividing the root edge gives one additional edge
        orientingNode = readyNodes.pop()
        orientingNodeInDeg = diNetwork.in_degree(orientingNode)
        if orientingNodeInDeg>2 or (orientingNodeInDeg==2 and orientingNode not in reticulations):
            network.add_edges_from([rootEdge])
            return False            
        children = set(network.neighbors(orientingNode))-set(diNetwork.predecessors(orientingNode))
        newArcs = []
        for child in children:
            newArcs +=[(orientingNode,child)]
            if diNetwork.in_degree(child)==1 or child not in reticulations:
                readyNodes.add(child)
            numberOrientedEdges+=1
        diNetwork.add_edges_from(newArcs)

    network.add_edges_from([rootEdge])
    if numberOrientedEdges < numberNetworkEdges+1:
        return False
    return diNetwork

def RootAtEdge(network,rootEdge,ClassChecker=ClassAllNetworks):
    noOfReticulations = len(network.edges)-len(network.nodes)+1
    nonLeafNodes = []
    for node in network.nodes:
        if network.degree(node)>1:
            nonLeafNodes +=[node]
    for reticulations in itertools.combinations(nonLeafNodes,noOfReticulations):
        result = OrientationAlgorithmBinary(network,rootEdge,set(reticulations))
        if result and ClassChecker(result):
            return reticulations
    return False      

def ClassRootableStupid(network,ClassChecker=ClassAllNetworks):
    rootings = dict()
    for e in network.edges:
        result = RootAtEdge(network,e,ClassChecker)
        if result:
            rootings[e]=result
    return rootings


###############################################
############   MORE ORIENTATION  ##############
###############################################

#Produces the generator includingDeg2Nodes
def GeneratorWithDeg2(network):
    generator = deepcopy(network)
    todoNodes = list(generator.nodes)
    while todoNodes:
        current = todoNodes.pop()
        if generator.degree(current)==1:
            generator.remove_node(current)
    return generator    

#Finds the side of the generator with deg2 nodes that starts with a given generator node and an adjacent node
def CompleteSide(generator,generatorNode,sideNode):
    previous = generatorNode
    current = sideNode
    side = (previous,current,)
    while generator.degree(current)==2:
        for nb in generator.neighbors(current):
            if nb != previous:
                next = nb
        previous,current = current,next
        side+=(current,)
    return side
        
#Returns all sides of the network
def FindSides(network):
    generator = GeneratorWithDeg2(network)    
    todoNodes = list(generator.nodes)
    sides=set()
    #Do something separate for cycles: in those cases, the sides do not end in generator nodes
    numberOfDeg3Nodes = 0
    while todoNodes:
        current = todoNodes.pop()
        if generator.degree(current)==3:
            numberOfDeg3Nodes +=1
            #Now look for sides around the generator node
            for nb in network.neighbors(current):
                side = CompleteSide(generator,current,nb)
                if not tuple(reversed(side)) in sides:
                    sides.add(side)
    if numberOfDeg3Nodes == 0:
        cycle = nx.cycle_basis(generator)[0]
        side = tuple([cycle[-1]]+cycle+[cycle[0]])
        sides.add(side)
    return sides
    
        
#Reduces chains in the network using the set of sides
#When reducing, it removes internal nodes of the chain.
def ReduceChains(network, length):
    sides = FindSides(network)
    sidesDict = dict()
    reducedNetwork = nx.Graph()
    for side in sides:
        if len(side)>length+2:
            #Keep the first l-1 and the last leaf on the chain
            newSide = side[:length]+side[-2:]
        else:
            newSide = side[:]
        sidesDict[newSide]=side
        #Now add the new side to the new network    
        previousNode = newSide[0]
        for node in newSide[1:]:
            reducedNetwork.add_edge(previousNode, node)
            previousNode = node
    #Add the leaves to the reduced network
    for node in list(reducedNetwork.nodes):
        if reducedNetwork.degree(node)==2:
            for nb in network.neighbors(node):
                if network.degree(nb)==1:
                    leaf = nb
            reducedNetwork.add_edge(node,leaf)
    return reducedNetwork, sidesDict

#Determines all root-edges of network
#Uses ClassChecker to check whether an orientation is in the desired class
#Uses that the class is length-chain reducible, blob-determined, and leaf-addable
def ReductionRooting(network,length,ClassChecker=ClassAllNetworks):
    redNw,sidesDict = ReduceChains(network,length)
    redRootings = ClassRootableStupid(redNw,ClassChecker)
    rootings = dict()
    for redSide,side in sidesDict.items():
        if len(side)<length+3:
        #Go through all leaf-edges and internal edges and copy the rootability to the network
            previous = redSide[0]
            for current in redSide[1:]:
                #Check if there is a leaf edge attached to current node
                for nb in redNw.neighbors(current):
                    nbDegree = redNw.degree(nb)
                    if nbDegree==1: 
                        #Check the leaf edge (current,nb)
                        rootingAtEdge = redRootings.get((current,nb))
                        if not rootingAtEdge:
                            rootingAtEdge = redRootings.get((nb,current))
                        if rootingAtEdge:
                            rootings[(current,nb)] = rootingAtEdge
                #Check the internal edge (previous,current)
                rootingAtEdge = redRootings.get((previous,current))
                if not rootingAtEdge:
                    rootingAtEdge = redRootings.get((current,previous))
                if rootingAtEdge:
                    rootings[(previous,current)] = rootingAtEdge
                previous = current            
        else:
        #Infer rootings from leaf edge rootability.
            n = len(side) - 2
            for index,current in enumerate(redSide[1:-1]):
                i=index+1
                #Find the leaf edge attached to current node
                for nb in redNw.neighbors(current):
                    nbDegree = redNw.degree(nb)
                    if nbDegree==1: 
                        leaf = nb
                #Check the leaf edge (current,leaf)
                rootingAtEdge = redRootings.get((current,leaf))
                if not rootingAtEdge:
                    rootingAtEdge = redRootings.get((leaf,current))
                if rootingAtEdge:
                    #If rootable at leaf edge, infer rootings of the original network
                    #First the leaf edges on the side
                    for j in range(i,n-(length-i)+1):
                        #Find the leaf edge at position j
                        sideNode = side[j]
                        for nb in network.neighbors(sideNode):
                            if network.degree(nb)==1:
                                rootings[(sideNode,nb)] = rootingAtEdge
                    #Now the internal edges of the side
                    for j in range(i-1,n-(length-i)+1):
                        rootings[(side[j],side[j+1])] = rootingAtEdge
    return rootings

        

#Determines one Class-rootedge of the network if it exists (False otherwise)
#Uses ClassChecker to check whether an orientation is in the desired class
#Uses that the class is length-chain reducible, blob-determined, and leaf-addable
def LevelStuff(network,length,ClassChecker=ClassAllNetworks):
    #Calculate the blobs
    blobs = list(nx.biconnected_component_subgraphs(network))
    #Prepare the partially directed network that we will condense into T_CN
    partiallyOrientedNetwork = network.to_directed()
    #Empty list to store all orientations for all blobs
    blobOrientations = []
    #For each blob, compute the orientations
    for blob in blobs:
        if len(blob)>2: 
        #Add leaves to degree 2 nodes in the biconnected component first, to make actual blobs               
            leafEdges = []
            for node in blob.nodes:
                if blob.degree(node)==2:
                    for nb in network.neighbors(node):
                        if nb not in blob:
                            leafEdges+=[(node,nb)]
            blob.add_edges_from(leafEdges)
        #Find all rootings of the blob
            rootingsBlob = ReductionRooting(blob,length,ClassChecker)
            blobOrientations+=[rootingsBlob]
            if not rootingsBlob:
                #If there is no orientation for this blob, then there is no orientation for the network.
                return False
        #Partially orient at the leaves, accoriding to where the blob can be rooted
            for leafEdge in leafEdges:
                if not (leafEdge in rootingsBlob or (leafEdge[1],leafEdge[0]) in rootingsBlob):
                    #If now both arcs are gone, there is no rooting of the original network, so we may return False
                    if not partiallyOrientedNetwork.has_edge(*leafEdge):
                        return False
                    partiallyOrientedNetwork.remove_edge(leafEdge[1],leafEdge[0])
        #For trivial biconnected components, add a trivial list to the blob orientations, so that indices still match between blobOrientations, and blobs
        else:
            blobOrientations+=[[]]
    #Create T_CN by condensing partiallyOrientedNetwork
    T_CN=nx.condensation(partiallyOrientedNetwork)
    #Find the root of T_CN if it exists; if it does not, the network is not C-orientable
    rootComponent = IsRootedTree(T_CN)
    if not type(rootComponent)==int:
        return False
        
    #Go through all edges to find all orientations    
    rootings = dict()
    rootComponentNodes = T_CN.nodes(data=True)[rootComponent]['members']
    for rootEdge in network.edges:
        reticulations = []
        edgesToContinueAt = False
        #Check if the edge is in the rootComponent
        if rootEdge[0] in rootComponentNodes and rootEdge[1] in rootComponentNodes:
            #Check if the edge is a root edge of one of the blobs
            for i,blob in enumerate(blobs):
                if blob.has_edge(*rootEdge):
                    if len(blob)==2:
                        edgesToContinueAt = set([rootEdge,(rootEdge[1],rootEdge[0])])
                        break
                    elif rootEdge in blobOrientations[i]: 
                        reticulations+=blobOrientations[i][rootEdge]
                        edgesToContinueAt = LeafEdges(blob)
                        break
                    elif (rootEdge[1],rootEdge[0]) in blobOrientations[i]:
                        reticulations+=blobOrientations[i][(rootEdge[1],rootEdge[0])]
                        edgesToContinueAt = LeafEdges(blob)
                        break
            #If it is a root edge, continue finding the whole orientation
            if edgesToContinueAt:     
                #Continue to root all other blobs, by moving away from the blob with the root.
                #edgesToContinueAt keeps a list of edges along which we still have to move away from the root
                while edgesToContinueAt:
                    edge = edgesToContinueAt.pop()
                    for i,blob in enumerate(blobs):
                        if edge[1] in blob :
                            #Continue at trivial biconnected components with an endpoint edge[1] (but not edge[0], to prevent cycling in the algorithm)
                            if len(blob)==2 and edge[0] not in blob:
                                otherNode = False
                                for v in blob:
                                    if v!= edge[1]:
                                        otherNode = v
                                edgesToContinueAt.add((edge[1],otherNode))
                            #Continue at blobs that contain edge[1] in the interior (so the degree of edge[1] in the blob is not 1)
                            elif len(blob)>2 and blob.degree(edge[1])!=1:
                                edgesToContinueAt|= LeafEdges(blob)-set([(edge[1],edge[0])])
                                if edge in blobOrientations[i]:
                                    reticulations+=blobOrientations[i][edge]
                                else:
                                    reticulations+=blobOrientations[i][(edge[1],edge[0])]
                                break
                rootings[rootEdge]=reticulations
    return rootings

            
            
            
            
                   

#Returns a list of leaf edges in the order (neighbor,leaf) for each edge
def LeafEdges(network):
    leafEdges=set()
    for node in network.nodes:
        if network.degree(node)==1:
            for nb in network.neighbors(node):
                leafEdges.add((nb,node))
    return leafEdges
                        
            




