import community
import networkx as nx
import matplotlib.pyplot as pl
from copy import deepcopy


class Status(object):
    def __init__(self):
        self.com2nodes = dict([])
        self.node2com = dict([])
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])
        self.n_edges = 0

    def __str__(self):
        return ("node2com:"
                + str(self.node2com)
                + '\ncom2nodes: '
                + str(self.com2nodes)
                + "\ndegrees: "
                + str(self.degrees)
                + "\ninternals: "
                + str(self.internals))

    def copy(self):
        '''
           Perform a deep copy of status.
        '''
        new_status = Status()
        new_status.node2com = self.node2com
        new_status.com2nodes = self.com2nodes
        new_status.internals = self.internals
        new_status.degrees = self.degrees
        new_status.gdegrees = self.gdegrees
        new_status.n_edges = self.n_edges
        new_status.com_mods = self.com_mods
        return new_status

    def init(self, graph, part=None):
        count = 0
        self.com2nodes = dict([])
        self.node2com = dict([])
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])
        self.n_edges = graph.size()
        self.com_mods = {}

        for node in graph.nodes():
            self.node2com[node] = count

            if not count in self.com2nodes:
                self.com2nodes[count] = set()
                self.com_mods[count] = 0

            self.com2nodes[count].add(node)

            deg = float(graph.degree(node))
            if deg < 0:
                error = "Bad graph type ({})".format(type(graph))
                raise ValueError(error)
            self.degrees[count] = deg
            self.gdegrees[node] = deg
            edge_data = graph.get_edge_data(node, node, {})
            self.loops[node] = 0
            self.internals[count] = 0
            count += 1

class State:
    def __init__(self):
        self.com_t = None
        self.com_s = None
        self.internals_t = None
        self.internals_s = None
        self.degrees_t = None
        self.degrees_s = None
        self.mod_s = None
        self.mod_t = None
        self.nodes_s = None
        self.nodes_t = None
        self.moved = []

    def reverse(self, status):
        status.internals[self.com_t] = self.internals_t
        status.internals[self.com_s] = self.internals_s
        status.degrees[self.com_t] = self.degrees_t
        status.degrees[self.com_s] = self.degrees_s
        status.com_mods[self.com_s] = self.mod_s
        status.com_mods[self.com_t] = self.mod_t
        status.com2nodes[self.com_t] = self.nodes_t
        status.com2nodes[self.com_s] = self.nodes_s

        for v in self.moved:
            status.node2com[v] = self.com_s
            
#@profile        
def merge_coms(l_comA, l_comB, G, status, state=None, verb=False):
    if verb:
       print('\n\nMerging: {} + {}'.format(l_comA, l_comB))
    
    if l_comA == l_comB:
        return l_comA

    comA = status.com2nodes[l_comA]
    comB = status.com2nodes[l_comB]

    if verb:
       print('Communities:\n1: {}\n2: {}'.format(comA, comB))

    if len(comA) > len(comB):
        target = comA
        t_label = l_comA
        source = comB
        s_label = l_comB
    else:
        target = comB
        t_label = l_comB
        source = comA
        s_label = l_comA

    if state:
        state.com_t = t_label
        state.com_s = s_label
        state.internals_s = status.internals[s_label]
        state.internals_t = status.internals[t_label]
        state.degrees_s = status.degrees[s_label]
        state.degrees_t = status.degrees[t_label]
        state.mod_s = status.com_mods[s_label]
        state.mod_t = status.com_mods[t_label]
        state.nodes_s = status.com2nodes[s_label].copy()
        state.nodes_t = status.com2nodes[t_label].copy()   

    diff = source - target

    for v in diff:
        degree = G.degree(v)
        in_degree = 0

        for n in G.neighbors(v):
            in_degree += 1 if status.node2com[n] == t_label else 0

        status.internals[t_label] += in_degree
        status.degrees[t_label] += degree
        status.node2com[v] = t_label
        status.com2nodes[t_label].add(v)

        if state:
            state.moved.append(v)

    del status.internals[s_label]
    del status.degrees[s_label]
    del status.com2nodes[s_label]
    del status.com_mods[s_label]

    return t_label


def cliq2com(cliq, G, status):
    coms = list(set([status.node2com[v] for v in cliq]))
    new_com = coms[0]

    for c in coms:
        new_com = merge_coms(new_com, c, G, status)

    return new_com

        
#@profile
def modularity(status, verbose=0, changed=None):
    changed = changed if changed else set(status.node2com.values())
    links = float(status.n_edges)
    result = 0.

    if verbose:
        print('\n     Calculating modularity...')
        print('     Number of edges in a graph: {}'.format(links))

#    for community in set(status.node2com.values()):
#        in_degree = status.internals.get(community, 0.)
#        degree = status.degrees.get(community, 0.)
#
#        if verbose > 1:
#            print('     Community: {}, in-degree: {}, degree: {}'.format(community, in_degree, degree))
#
#        if links > 0:
#            result += in_degree / links - ((degree / (2. * links)) ** 2)
        
    for community in changed:
        in_degree = status.internals.get(community, 0.)
        degree = status.degrees.get(community, 0.)
        com_mod = 0

        if verbose > 1:
            print('     Community: {}, in-degree: {}, degree: {}'.format(community, in_degree, degree))

        if links > 0:
            com_mod += in_degree / links - ((degree / (2. * links)) ** 2)

        status.com_mods[community] = com_mod

#    if len(status.com_mods.keys()) != len(status.com2nodes.keys()):
#        raise Exception('Mismatch: {} /= {}'.format(status.com_mods.keys() , status.com2nodes.keys()))

    if verbose:
        print('     Result: {}\n'.format(result))

    modularity = sum(status.com_mods.values())

#    if abs(modularity - result) > 0.000001:
#        print("ERROR: {} /= {}".format(modularity, result))

    return modularity


#@profile
def cluster(G):
    max_mod = -1
    touched = set()
    status = Status()
    status.init(G)
    _status = status.copy()

    polys = [frozenset(c) for c in nx.find_cliques(G) if len(c) > 2]
    polys.sort(key=lambda p: -len(p))

    print('Phase 1... {} cliques'.format(len(polys)))

    for p in polys:
        inter = 0

        for v in p:
            inter += 1 if v in touched else 0

        if len(p) - inter + 1 < 3:
            continue

        touched = touched.union(p)

        new_com = cliq2com(p, G, _status)
        new_mod = modularity(_status, 0, [new_com])

        if new_mod > max_mod:
            max_mod = new_mod
            status = _status.copy()
        else:
            continue

    min_del = 0.0000001
    merged = set()
    used = set()
    coms = set(status.node2com.values())
    count = 0

    print('Phase 2... {} communities'.format(len(coms)))

    for com in coms:
        # Check if community exists.
        if com in merged:
            continue

        for com_ in coms:
            if com == com_ or com_ in merged or frozenset([com, com_]) in used:
                continue

            used.add(frozenset([com, com_]))
            count += 1

            #print('Merging: {} + {}'.format(com, com_))

            _status = status
            state = State()
            #print('\nBefore merge: {}'.format(_status))
            new_com = merge_coms(com, com_, G, _status, state, verb=False)
            new_mod = modularity(_status, 0, [new_com])

            if new_mod > max_mod:
                max_mod = new_mod
                status = _status
                merged.add(com + com_ - new_com)

                if min_del > (new_mod - max_mod):
                    print('Merged {} times.'.format(count))
                    print('Final modularity: {}'.format(modularity(status, 0)))
                    # print('Final partition: {}'.format(relabel(status.node2com)))
                    return (status.node2com, G)
                
                break
            else:
                print('Merged {} times.'.format(count))
                print('Final modularity: {}'.format(modularity(status, 0)))
                # print('Final partition: {}'.format(relabel(status.node2com)))
                return (status.node2com, G)
                #print("EHHHH...")
                #print('After merge: {}'.format(_status))
                #print('Reversing... ')
                state.reverse(_status)
                #print('After reverse: {}'.format(_status))
                
    print('Merged {} times.'.format(count))
    print('Final modularity: {}'.format(modularity(status, 0)))
    # print('Final partition: {}'.format(relabel(status.node2com)))
    return (status.node2com, G)
