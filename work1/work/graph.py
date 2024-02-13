class DirectedGraph:
    """Abstract base class for a directed graph.

    A functional directed graph class can be obtained by inheriting from 
    this class and overriding the methods has_edge and add_edge.  All other
    methods have default implementations, which may not be the most efficient.
    These other methods should also be overriden as appropriate to improve
    efficiency.
    """
    def __init__(self, num_vertices):
        """Constructs a directed graph with num_vertices vertices and zero edges"""
        self._num_vertices = num_vertices
    
    def has_edge(self, i, j):
        """Returns True if the graph contains the directed edge (i, j), False otherwise."""
        raise NotImplementedError
        
    def add_edge(self, i, j):
        """Adds the directed edge (i, j) to the graph."""
        raise NotImplementedError
        
    def out_edges(self, i):
        """Returns a list of directed edges outgoing from vertex i."""
        return [(i, j) for j in range(self._num_vertices) if self.has_edge(i, j)]
    
    def in_edges(self, j):
        """Returns a list of directed edges incoming to vertex j."""
        return [(i, j) for i in range(self._num_vertices) if self.has_edge(i, j)]
    
    def outdegree(self, i):
        """Returns the outdegree of vertex i."""
        return len(self.out_edges(i))
    
    def indegree(self, i):
        """Returns the indegree of vertex i."""
        return len(self.in_edges(i))
    
    def degree(self, i):
        """Returns the degree of vertex i."""
        return self.indegree(i) + self.outdegree(i)
        
    def add_edges(self, edges):
        """Adds all edges from a list to the graph."""
        for i, j in edges:
            self.add_edge(i, j)
            
    def num_vertices(self):
        """Returns the number of vertices in the graph."""
        return self._num_vertices

    def num_edges(self):
        """Returns the number of edges in the graph."""
        return len(tuple(self.edges()))
    
    def edges(self):
        """Returns an iterator over the edges of the graph."""
        for i in range(self._num_vertices):
            for edge in self.out_edges(i):
                yield edge
    
    def __str__(self):
        """Returns a string representation of the graph, so that it may be printed."""
        return "DirectedGraph with %d vertices and %d edge(s):\n%s" % (self.num_vertices(),
                                                                       self.num_edges(),
                                                                       sorted(self.edges()))

class TrivialSetDirectedGraph(DirectedGraph):
    """A trivial implementation of a Directed Graph that simply stores edges in a set.
    Not meant for serious use.
    """
    
    def __init__(self, num_vertices):
        # call the next constructor in Python's Method Resolution Order
        super().__init__(num_vertices)
        # start with an empty set of edges
        self._edges = set()
    
    def has_edge(self, i, j):
        return (i, j) in self._edges
        
    def add_edge(self, i, j):
        self._edges.add((i, j))
        
class AdjacencyListDirectedGraph(DirectedGraph):
    """An implementation of a Directed Graph that uses adjacency lists to store edges."""
        
    def __init__(self, num_vertices):
        super().__init__(num_vertices)
        self._out_lists = [[] for i in range(num_vertices)]
        self._in_lists = [[] for i in range(num_vertices)]
    
    def add_edge(self, i, j):
        self._out_lists[i].append(j)
        self._in_lists[j].append(i)
    
    def has_edge(self, i, j):
        return j in self._out_lists[i]
        
    def out_edges(self, i):
        return [(i, j) for j in self._out_lists[i]]
        
    def in_edges(self, j):
        return [(i, j) for i in self._in_lists[j]]
    
    def indegree(self, i):
        return len(self._in_lists[i])
        
    def outdegree(self, i):
        return len(self._out_lists[i])
    
    def cycle(self, i ,j):
        if len(self._out_lists[j]) == 0:
            return False
        outver = self._out_lists[j][0]
        while(1):
            if outver == i:
                return True
            if len(self._out_lists[outver]) == 0:
                break
            outver = self._out_lists[outver][0]
        return False
            
    
    

class VertexLabeledDirectedGraph(DirectedGraph):
    """A mixin class that allows for vertices in a directed graph to have labels."""
    
    def __init__(self, num_vertices):
        # call the next constructor in Python's Method Resolution Order
        super().__init__(num_vertices)
        self._vertex_labels = [None] * num_vertices

    def set_vertex_label(self, i, label):
        """Adds a label to vertex i."""
        self._vertex_labels[i] = label
       
    def vertex_label(self, i):
        """Returns the label of vertex i or None if it is not assigned a label"""
        return self._vertex_labels[i]
    
class VertexLabeledALDirectedGraph(AdjacencyListDirectedGraph, VertexLabeledDirectedGraph):
    pass

def find_cycle_in_eulerian_graph(g, start_vertex):
    """Finds an arbitrary cycle starting at a given vertex in a directed graph.
    Assumes that the graph is Eulerian, and thus that such a cycle exists."
    
    Args:
        g: a DirectedGraph object
        start_vertex: the index of the vertex from which the cycle should start
    Returns:
        A cycle represented by a list of indices of the vertices traversed, in order, by the cycle.  
        The first and last entries of the list will be identical and equal to start_vertex.
    """
    vertex = start_vertex
    ans = [start_vertex]
    s = set()
    delList = []
    for i in range(g.num_edges()):
        verList = set(g.out_edges(vertex))
        for j in verList:
            if j in delList:
                continue
            elif j[1] == vertex:
                continue
            else:
                ans.append(j[1])
                vertex = j[1] 
                delList.append(j)
                break
    return ans

def is_cycle(g, path, vertex):
    """Returns True if path is a cycle in g that starts and ends with vertex."""
    # check for cyclical path with given start/end vertex
    if len(path) < 2 or path[0] != vertex or path[-1] != vertex:
        return False
    # check for distinct edges
    path_edges = set(zip(path[:-1], path[1:]))    
    if len(path_edges) != len(path) - 1:
        return False
    # check for existence of edges
    for i, j in path_edges:
        if not g.has_edge(i, j):
            return False
    return True