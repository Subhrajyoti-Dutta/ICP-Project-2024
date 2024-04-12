import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
// This class represents a node in a graph data structure
class Node {

    // The position of the node
    private final int pos;

    // List of edges connected to this node
    private final List<Edge> edges;

    // List of neighboring nodes connected by the edges
    private final List<Node> neighbours;

    // Constructor to initialize the node with its position
    public Node(int pos) {
        this.pos = pos;
        this.edges = new ArrayList<>();  // Initialize an empty list for edges
        this.neighbours = new ArrayList<>(); // Initialize an empty list for neighbours
    }

    // Getter method to access the position of the node
    public int getPos() {
        return pos;
    }

    // Method to add a new edge and its corresponding neighbour node
    public void addNeighbour(Edge edge, Node node) {
        edges.add(edge);  // Add the edge to the list
        neighbours.add(node);  // Add the neighbour node to the list
    }

    // Getter method to retrieve the list of edges connected to this node
    public List<Edge> getEdges() {
        return edges;
    }

    // Getter method to retrieve the list of neighbour nodes connected to this node
    public List<Node> getNeighbours() {
        return neighbours;
    }

    // Override the toString method to provide a human-readable representation of the node
    @Override
    public String toString() {
        return "Node(" + pos + ')';
    }
}
// This class represents an edge in a graph data structure
class Edge {

    // Two nodes that this edge connects
    private final Node[] nodes;

    // Weight associated with the edge (cost, distance, etc.)
    private final int weight;

    // Constructor to initialize the edge with two connected nodes and weight
    public Edge(Node node1, Node node2, int weight) {
        this.nodes = new Node[]{node1, node2};
        this.weight = weight;
    }

    // Getter method to access the weight of the edge
    public int getWeight() {
        return weight;
    }

    // Method to get the neighbouring node of a given node relative to this edge
    public Node getNeighbour(Node node) {
        return (node == nodes[0]) ? nodes[1] : nodes[0];
    }

    // Getter method to access a specific node based on its index (0 or 1) in the edge
    public Node getNode(int index) {
        return nodes[index];
    }

    // Override the toString method to provide a human-readable representation of the edge
    @Override
    public String toString() {
        return nodes[0] + "--" + weight + "--" + nodes[1];  // Show connected nodes and weight
    }
}

// This class represents a graph data structure
class Graph {

    // List of nodes in the graph
    protected final List<Node> nodes;

    // List of edges connecting the nodes
    protected final List<Edge> edges;

    // Number of nodes in the graph
    protected final int size;

    // Parent array used for efficient union-find operations (for minimum spanning tree algorithms)
    protected final int[] parent;

    // Constructor to initialize the graph with a specified size
    public Graph(int size) {
        this.size = size;
        this.edges = new ArrayList<>();
        this.nodes = new ArrayList<>();
        this.parent = new int[size];

        // Initialize nodes and parent array (each node initially belongs to its own component)
        for (int i = 0; i < size; i++) {
            nodes.add(new Node(i));
            parent[i] = i;
        }
    }

    // Getter method to access the number of nodes in the graph
    public int getSize() {
        return size;
    }

    // Getter method to access the number of edges in the graph
    public List<Edge> getEdges() {
        return edges;
    }

    // Method to add a new edge between two nodes with a weight
    public void addEdge(int node1, int node2, int weight) {
        if (node1 == node2) {
            return;  // Don't allow self-loops
        }

        // Get references to the nodes
        Node getNode1 = nodes.get(node1);
        Node getNode2 = nodes.get(node2);

        // Create a new edge object
        Edge edge = new Edge(getNode1, getNode2, weight);

        // Add the edge to the graph and update the neighbor lists of the nodes
        edges.add(edge);
        getNode1.addNeighbour(edge, getNode2);
        getNode2.addNeighbour(edge, getNode1);

        // Perform union-find operation to merge connected components (relevant for minimum spanning tree algorithms)
        union(node1, node2);
    }

    // Helper method to find the root of a node in the union-find data structure (path compression)
    protected int find(int node) {
        if (parent[node] != node) {
            parent[node] = find(parent[node]);  // Recursively find the root and compress the path
        }
        return parent[node];
    }

    // Helper method to merge two connected components using union-find
    protected void union(int node1, int node2) {
        int root1 = find(node1);
        int root2 = find(node2);
        if (root1 != root2) {
            parent[root1] = root2;  // Attach the smaller component's root to the larger component's root
        }
    }

    // Method to check if all nodes in the graph are connected (useful for finding minimum spanning trees)
    public boolean isConnected() {
        return nodes.stream()
                .allMatch(node -> find(node.getPos()) == find(nodes.getFirst().getPos()));
    }

    // Override the toString method to provide a human-readable representation of the graph (shows edges)
    @Override
    public String toString() {
        return edges.stream()
                .map(Objects::toString)
                .collect(Collectors
                        .joining(System.lineSeparator()));
    }

    public List<Node> getNodes() {
        return nodes;
    }

}

class AcyclicGraph extends Graph {

    public AcyclicGraph(int size) {
        super(size);
    }

    @Override
    public void addEdge(int node1, int node2, int weight) {
        if (!checkCyclic(node1, node2)) {
            super.addEdge(node1, node2, weight);
        }
    }

    private boolean checkCyclic(int start, int end) {
        Set<Node> visited = new HashSet<>();
        Node startNode = nodes.get(start);
        Node endNode = nodes.get(end);
        return hasCycle(startNode, endNode, visited);
    }

    private boolean hasCycle(Node current, Node target, Set<Node> visited) {
        if (current == target) {
            return true;
        }
        visited.add(current);
        return current.getNeighbours()
                .stream()
                .anyMatch(neighbour -> !visited
                        .contains(neighbour)
                        && hasCycle(neighbour, target, visited));
    }
}

class KruskalGraph extends AcyclicGraph {
    public KruskalGraph(Graph graph) {
        super(graph.getSize());

        PriorityQueue<Edge> minHeap = new PriorityQueue<>(Comparator.comparingInt(Edge::getWeight));
        minHeap.addAll(graph.edges);

        for (Edge edge : minHeap) {
            int root1 = find(edge.getNode(0).getPos());
            int root2 = find(edge.getNode(1).getPos());

            if (root1 != root2) {
                super.addEdge(edge.getNode(0).getPos(), edge.getNode(1).getPos(), edge.getWeight());
                union(root1, root2);
            }
        }
    }
}

class PrimGraph extends AcyclicGraph {

    public PrimGraph(Graph graph) {
        super(graph.getSize());
        Set<Node> unvisitedNodes = new HashSet<>(graph.getNodes());

        // Use priority queue to select the minimum weight edge efficiently
        PriorityQueue<Edge> minHeap = new PriorityQueue<>(Comparator.comparingInt(Edge::getWeight));

        // Choose a random starting node
        Node startNode = unvisitedNodes.iterator().next();
        unvisitedNodes.remove(startNode);

        // Add all edges of the start node to the priority queue
        minHeap.addAll(startNode.getEdges());

        while (!minHeap.isEmpty()) {
            Edge minEdge = minHeap.poll();
            Node current = minEdge.getNode(0);
            Node neighbour = minEdge.getNode(1);

            // If neighbour is unvisited, add the edge and update
            if (unvisitedNodes.remove(neighbour)) {
                super.addEdge(current.getPos(), neighbour.getPos(), minEdge.getWeight());
                minHeap.addAll(neighbour.getEdges().stream()
                        .filter(edge -> unvisitedNodes.contains(edge.getNeighbour(neighbour)))
                        .toList());
            }
        }
    }
}

public class Kruskal {
    public static void main(String[] args) throws IOException {
        int upper = 5100;
        Random random = new Random();
        Graph myGraph;
        File file = new File("time.csv");
        FileWriter fr = new FileWriter(file, false);
        File file2 = new File("time2.csv");
        FileWriter fr2 = new FileWriter(file2, false);
        int node1, node2;
        int rep = 3;
        for (int n = 100; n < upper; n += 100) {
            long startTime2 = 0, endTime2 = 0;
            long startTime = 0, endTime = 0;
            System.out.println("Graph Size = " + n);
            for (int j = 0; j < rep; j++) {
                myGraph = new Graph(n);
                int i = 0;
                do {
                    node1 = random.nextInt(n);
                    node2 = random.nextInt(n);
                    if (node1 == node2)
                        continue;
                    i++;
                    myGraph.addEdge(node1, node2, random.nextInt(100));
                } while (i < n * (n - 1) / 4);
                do {
                    node1 = random.nextInt(n);
                    node2 = random.nextInt(n);
                    if (node1 == node2)
                        continue;
                    myGraph.addEdge(node1, node2, random.nextInt(100));

                } while (!myGraph.isConnected());

                startTime += System.nanoTime();
                new PrimGraph(myGraph);
                endTime += System.nanoTime();

                startTime2 += System.nanoTime();
                new KruskalGraph(myGraph);
                endTime2 += System.nanoTime();
            }
            System.out.println("Total-time  = " + (endTime - startTime) / 1e9 / rep);
            fr.write(n + "," + (endTime - startTime) / 1e9 / rep + System.lineSeparator());

            System.out.println("Total-time  = " + (endTime2 - startTime2) / 1e9 / rep);
            fr2.write(n + "," + (endTime2 - startTime2) / 1e9 / rep + System.lineSeparator());
        }
        fr.close();
        fr2.close();
    }
}
