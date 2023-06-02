from rooting import *
import networkx as nx
import os
import sys
import csv

errorInFile = False

#####################################################################
#####################################################################
##############          Terminal arguments                ###########
#####################################################################
#####################################################################

#Set default values
option_out = False
option_out_argument = ""
option_file = False
option_file_argument = ""
option_class = False
option_class_argument = ""
option_simple = False
option_help = False

#Read the arguments
i = 1
while i < len(sys.argv):
    arg= sys.argv[i]
    if arg == "-f" or arg == "--file":
        option_file = True
        i+=1
        option_file_argument = sys.argv[i]
    if arg == "-o" or arg == "--output":
        option_out = True
        i+=1
        option_out_argument = sys.argv[i]
    if arg == "-c" or arg == "--class":
        option_class = True
        i+=1
        option_class_argument = sys.argv[i]
    if arg == "-s" or arg == "--simple":
        option_simple = True
    if arg == "-h" or arg == "--help":
        option_help = True
    i += 1



#Set the class parameters corresponding to the option_class_argument
nwClass = option_class_argument
ClassChecker = ClassAllNetworks
length = 1
if nwClass == "TC":
    ClassChecker = ClassTreeChild
    length = 3
elif nwClass == "SF":
    ClassChecker = ClassStackFree
    length = 3
elif nwClass == "O":
    ClassChecker = ClassOrchard
    length = 3
elif nwClass == "TB":
    ClassChecker = ClassTreeBased
    length = 2
elif not option_class:
    print("Error: have to choose a class!\n\n")
    option_help = True    
else:
    print("Error: invalid class chosen!\n\n")
    option_help = True



#Output the help text to the terminal if no argument is given, or if the help option is chosen.
if len(sys.argv)==1 or option_help:
    print("Mandatory arguments:\n -f or --file followed by the input file with a list of edges. One edge per line, vertices are positive integers separated by a comma\n -c or -class followed by a class:\n   TC: tree-child\n   TB: tree-based\n   SF: stack-free\n   O: Orchard \n\nOptional arguments:\n -o or --output followed by the output file name\n -s or --simple for output only consisting of the root-edge and the reticulation nodes.")
    sys.exit()






#####################################################################
#####################################################################
##############             Read the input                 ###########
#####################################################################
#####################################################################




#Read each line of the input file with name set by "option_file_argument"
edges = []
f = open("./"+option_file_argument, "rt")
reader = csv.reader(f, delimiter=',', quotechar='|')
for row in reader:
    if len(row)==2:
        edges.append((int(row[0]),int(row[1])))
f.close()

#Build the unrooted network to be oriented
network = nx.Graph()
network.add_edges_from(edges)



#####################################################################
#####################################################################
##############            Find an orientation             ###########
#####################################################################
#####################################################################

orientations = LevelStuff(network,length,ClassChecker)

output = ""
if orientations:
    rootableEdges = orientations.keys()
    #print("rootable at: ", rootableEdges)
    #print("not at: ", set(network.edges)-set(rootableEdges))
    if not option_simple:
        output += "The root edges, reticulations and orientations are:"
    else:
        output += "The orientations of the network consist of the following root edges and reticulations:\r\n"
    for rootEdge, reticulations in orientations.items():
        #compute the actual orientation, instead of only the retculations.
        dinetwork = OrientationAlgorithmBinary(network,rootEdge,reticulations)
        if not option_simple:
            output += "\r\n\r\nroot edge:\r\n   "+str(rootEdge)+"\r\nreticulations:\r\n   "+str(reticulations)+"\r\norientation:\r\n   "+str(dinetwork.edges)
        else:
            output += "\r\n"+str(rootEdge)+" "+str(reticulations)
    if not option_simple:
        output+="\r\n\r\nThe network cannot be rooted at the edges:"
        for non_root in set(network.edges)-set(rootableEdges):
            output+="\r\n"+str(non_root)
else:
    output = "There is no orientation of the given network in the desired class."


#####################################################################
#####################################################################
##############                   Output                   ###########
#####################################################################
#####################################################################


#Open the output file for writing
if option_out:
    f= open(option_out_argument,"w+")
    f.write("The input network was:\r\n")
    f.write("  "+str(network.edges)+"\r\n")
    f.write("The class is:\r\n")
    f.write("  "+option_class_argument+"\r\n")    
    f.write(output)
    f.close()
else:
    print(output)

