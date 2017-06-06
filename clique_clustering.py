import community
import networkx as nx
import matplotlib.pyplot as plt


class Status(object):
    node2com = {}
    total_weight = 0
    internals = {}
    degrees = {}
    gdegrees = {}

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
        new_status.node2com = self.node2com.copy()
        new_status.com2nodes = self.com2nodes.copy()
        new_status.internals = self.internals.copy()
        new_status.degrees = self.degrees.copy()
        new_status.gdegrees = self.gdegrees.copy()
        new_status.n_edges = self.n_edges
        return new_status

    def init(self, graph, part=None):
        """Initialize the status of a graph with every node in one community"""
        count = 0
        self.com2nodes = dict([])
        self.node2com = dict([])
        self.degrees = dict([])
        self.gdegrees = dict([])
        self.internals = dict([])
        self.loops = dict([])
        self.n_edges = graph.size()

        if part is None:
            for node in graph.nodes():
                self.node2com[node] = count

                if not count in self.com2nodes:
                    self.com2nodes[count] = set()
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
        else:
            for node in graph.nodes():
                com = part[node]
                self.node2com[node] = com
                deg = float(graph.degree(node))
                self.degrees[com] = self.degrees.get(com, 0) + deg
                self.gdegrees[node] = deg
                inc = 0.
                for neighbor, datas in graph[node].items():
                    if part[neighbor] == com:
                        if neighbor == node:
                            inc += 1
                        else:
                            inc += 1 / 2.
            self.internals[com] = self.internals.get(com, 0) + inc


def modularity(status, verbose=0):
    links = float(status.n_edges)
    result = 0.

    if verbose:
        print('\n     Calculating modularity...')
        print('     Number of edges in a graph: {}'.format(links))

    for community in set(status.node2com.values()):
        in_degree = status.internals.get(community, 0.)
        degree = status.degrees.get(community, 0.)

        if verbose > 1:
            print('     Community: {}, in-degree: {}, degree: {}'.format(community, in_degree, degree))

        if links > 0:
            result += in_degree / links - ((degree / (2. * links)) ** 2)

    if verbose:
        print('     Result: {}\n'.format(result))

    return result


# @profile
def merge_coms(l_comA, l_comB, G, status):
    '''
        Merge two communites together and update Status.
    '''
    if l_comA == l_comB:
        return l_comA

    comA = status.com2nodes[l_comA]
    comB = status.com2nodes[l_comB]

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

    del status.internals[s_label]
    del status.degrees[s_label]
    del status.com2nodes[s_label]

    return t_label


# @profile
def cliq2com(cliq, G, status):
    coms = list(set([status.node2com[v] for v in cliq]))
    res = coms[0]

    for c in coms:
        res = merge_coms(res, c, G, status)

    return res


#@profile
def cluster(G):
    max_mod = -1
    touched = set()
    status = Status()
    status.init(G)
    _status = status.copy()
    polys = sorted([frozenset(c) for c in nx.find_cliques(G) if len(c) > 2], key=lambda p: -len(p))

    print('Phase 1...')

    for p in polys:
        inter = 0

        for v in p:
            inter += 1 if v in touched else 0

        if len(p) - inter + 1 < 3:
            continue

        touched = touched.union(p)

        cliq2com(p, G, _status)
        new_mod = modularity(_status, 0)

        if new_mod > max_mod:
            max_mod = new_mod
            status = _status.copy()
        else:
            continue

    merged = set()
    coms = set(status.node2com.values())

    print('Phase 2...')

    for com in coms:
        # Check if community exists.
        if com in merged:
            continue

        for com_ in coms:
            _status = status.copy()

            if com == com_ or com_ in merged:
                continue

            new_com = merge_coms(com, com_, G, _status)
            new_mod = modularity(_status, 0)

            if new_mod > max_mod:
                max_mod = new_mod
                status = _status
                merged.add(com + com_ - new_com)
                break

    print('Final modularity: {}'.format(modularity(status, 0)))
    # print('Final partition: {}'.format(relabel(status.node2com)))
    return (status.node2com, G)
