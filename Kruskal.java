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
        super(size);  // Call the constructor of the parent class Graph
    }

    @Override
    public void addEdge(int node1, int node2, int weight) {
        // Check if adding the edge would create a cycle before adding it
        if (!checkCyclic(node1, node2)) {
            super.addEdge(node1, node2, weight);  // Call the addEdge method from Graph
        }
    }

    // Check if there's a cycle starting from node1 and reaching node2
    private boolean checkCyclic(int start, int end) {
        Set<Node> visited = new HashSet<>();  // Keep track of visited nodes to detect cycles
        Node startNode = nodes.get(start);
        Node endNode = nodes.get(end);
        return hasCycle(startNode, endNode, visited);
    }

    // Recursive helper function to perform the actual cycle detection
    private boolean hasCycle(Node current, Node target, Set<Node> visited) {
        // If the current node is the target node, a cycle is found
        if (current == target) {
            return true;
        }

        // Mark the current node as visited
        visited.add(current);

        // Recursively check the neighbours of the current node
        // If any neighbour (not already visited) leads to the target node (or a cycle), return true
        return current.getNeighbours()
                .stream()
                .anyMatch(neighbour -> !visited.contains(neighbour) && hasCycle(neighbour, target, visited));
    }
}
class KruskalGraph extends AcyclicGraph {

    public KruskalGraph(Graph graph) {
        super(graph.getSize());  // Inherit the size from the provided graph

        // Priority queue to efficiently select edges with minimum weight
        PriorityQueue<Edge> minHeap = new PriorityQueue<>(Comparator.comparingInt(Edge::getWeight));
        minHeap.addAll(graph.getEdges());  // Add all edges from the provided graph to the queue

        // Kruskal's algorithm for finding Minimum Spanning Tree
        for (Edge edge : minHeap) {
            // Get the root nodes (represent connected components) of the two nodes in the edge
            int root1 = find(edge.getNode(0).getPos());
            int root2 = find(edge.getNode(1).getPos());

            // If the nodes belong to different connected components (not yet merged)
            if (root1 != root2) {
                // Add the edge to the MST (assuming super.addEdge ensures acyclic property)
                super.addEdge(edge.getNode(0).getPos(), edge.getNode(1).getPos(), edge.getWeight());

                // Merge the connected components using union-find (assuming union and find methods are implemented)
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
        // Maximum number of nodes in the graph (adjustable)
        final int upperLimit = 5000;

        // Random number generator for creating graph edges
        Random random = new Random();

        // Files for storing execution time measurements
        File primTimeFile = new File("time.csv");
        File kruskalTimeFile = new File("time2.csv");

        // Writers for the CSV files
        FileWriter primFileWriter = new FileWriter(primTimeFile, false); // Overwrite existing content
        FileWriter kruskalFileWriter = new FileWriter(kruskalTimeFile, false);

        int node1, node2;
        // Number of repetitions for each graph size to get average time
        int repetitions = 3;

        for (int numNodes = 100; numNodes <= upperLimit; numNodes += 100) {
            long startTimePrim = 0, endTimePrim = 0;
            long startTimeKruskal = 0, endTimeKruskal = 0;

            System.out.println("Graph Size = " + numNodes);

            for (int j = 0; j < repetitions; j++) {
                // Create a new graph with the specified number of nodes
                Graph myGraph = new Graph(numNodes);

                // Generate random edges until the graph is connected (ensures MST is possible)
                int edgesAdded = 0;
                do {
                    node1 = random.nextInt(numNodes);
                    node2 = random.nextInt(numNodes);
                    if (node1 == node2) {
                        continue;  // Skip self-loops
                    }
                    myGraph.addEdge(node1, node2, random.nextInt(100));
                } while (++edgesAdded < numNodes * (numNodes - 1) / 4 && !myGraph.isConnected());

                // Measure execution time for Prim's algorithm
                startTimePrim += System.nanoTime();
                new PrimGraph(myGraph);
                endTimePrim += System.nanoTime();

                // Measure execution time for Kruskal's algorithm
                startTimeKruskal += System.nanoTime();
                new KruskalGraph(myGraph);
                endTimeKruskal += System.nanoTime();
            }

            // Calculate and print average execution time for Prim's algorithm
            double avgPrimTime = (endTimePrim - startTimePrim) / 1e9 / repetitions;
            System.out.println("Average Time (Prim's) = " + avgPrimTime + " seconds");
            primFileWriter.write(numNodes + "," + avgPrimTime + System.lineSeparator());

            // Calculate and print average execution time for Kruskal's algorithm
            double avgKruskalTime = (endTimeKruskal - startTimeKruskal) / 1e9 / repetitions;
            System.out.println("Average Time (Kruskal's) = " + avgKruskalTime + " seconds");
            kruskalFileWriter.write(numNodes + "," + avgKruskalTime + System.lineSeparator());
        }

        // Close the file writers
        primFileWriter.close();
        kruskalFileWriter.close();
    }
}