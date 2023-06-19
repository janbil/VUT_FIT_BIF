from Bio import Phylo
from Bio import SeqIO
import csv

#Read ancestrals.csv and convert posterior probabilities to floats.
with open('ancestrals.csv', newline='') as csvfile:
    data = list(csv.reader(csvfile))
ii = 0
for i, row in enumerate(data):
    if(ii!=0):
        for j, val in enumerate(row):
            if val == '-':
                data[i][j] = 0
            else:
                data[i][j] = float(val)
    ii = i+1

#For given node picks the most probable aminoacid at each position.
#Input is number of the node and output is dictionary {position:aminoacid}
def lookup_ancestrals(node):
    res = {}
    for row in data:
        curr_max = 0
        curr_index = 0
        if row[0] == node:
            for i in range(2,len(row)):
                if row[i]>curr_max:
                    curr_max = row[i]
                    curr_index = i
            res[row[1]] = data[0][curr_index]
    return res

#For given node finds leaf nodes and distance to them
#Input Clade node and 0. Output List of tuples in form (name_of_leaf_node, depth)
def get_leaves(node, depth):
    res = []
    if len(node.clades) == 0:
        return [(node.name,depth)]

    for child in node.clades:
        new_depth = depth + child.branch_length
        for leaf in get_leaves(child, new_depth):
            res.append(leaf)
    return res

#Converts dictionary {pos:aminoacid} to string
def dict_to_str(dict):
    str = ""
    for pos in dict:
        str+=dict[pos]
    return str

#Read tree and msa.
tree = Phylo.read("tree.tre", "newick")
record_dict = SeqIO.to_dict(SeqIO.parse("msa.fasta", "fasta"))

# BFS and for each node and position computes sum of depths to all leaf nodes separetly for '-' and other letters
# If sum of depths at position where '-' occures is greater then sum of depths at position where '-' doesnt occure, change to '-'
queue = [tree.clade]
nodes = {}
while len(queue) > 0:
    cur_node = queue.pop(0)
    anc = lookup_ancestrals(cur_node.confidence)
    sequences = get_leaves(cur_node, 0)
    for pos in anc:
        space = 0
        nonspace = 0
        for sequence in sequences:
            letter = (record_dict[sequence[0]].seq._data).decode('UTF-8')[int(pos)-1]
            if letter == '-':
                space += sequence[1]
            else:
                nonspace += sequence[1]
        if(space> nonspace):
            #print('change ' + str(cur_node.confidence) + ' at pos ' + str(pos) + ' ,space = ' + str(space) + ' nonspace = ' + str(nonspace))
            anc[pos] = '-'
    nodes[cur_node.confidence] = anc
    if len(cur_node.clades)==2:
        if len(cur_node.clades[0].clades) > 0:
            queue.append(cur_node.clades[0])
        if len(cur_node.clades[1].clades) > 0:
            queue.append(cur_node.clades[1])
            
#Save sequences to right files
for num in nodes:
    name = "out/node_"+ str(num) + ".fas"
    file = open(name, "w")
    a = file.write(dict_to_str(nodes[num]))
    file.close()