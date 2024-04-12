//@formatter:on

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

class Node {
    private final int pos;
    private final List<Edge> edges;
    private final List<Node> neighbours;

    public Node(int pos) {
        this.pos = pos;
        this.edges = new ArrayList<>();
        this.neighbours = new ArrayList<>();
    }

    public int getPos() {
        return pos;
    }

    public void addNeighbour(Edge edge, Node node) {
        edges.add(edge);
        neighbours.add(node);
    }

    public List<Edge> getEdges() {
        return edges;
    }

    public List<Node> getNeighbours() {
        return neighbours;
    }

    @Override
    public String toString() {
        return "Node(" + pos + ')';
    }
}

class Edge {
    private final Node[] nodes;
    private final int weight;

    public Edge(Node node1, Node node2, int weight) {
        this.nodes = new Node[]{node1, node2};
        this.weight = weight;
    }

    public int getWeight() {
        return weight;
    }

    public Node getNeighbour(Node node) {
        if (node == nodes[0])
            return nodes[1];
        else
            return nodes[0];
    }


    public Node getNode(int index) {
        return nodes[index];
    }

    @Override
    public String toString() {
        return nodes[0] + "--" + weight + "--" + nodes[1];
    }
}

class Graph {
    protected final List<Node> nodes;
    protected final List<Edge> edges;
    protected final int size;
    protected final int[] parent; // Union-Find parent array
    public long time = 0;

    public Graph(int size) {
        this.size = size;
        this.edges = new ArrayList<>();
        this.nodes = new ArrayList<>();
        this.parent = new int[size];
        for (int i = 0; i < size; i++) {
            nodes.add(new Node(i));
            parent[i] = i; // Each node initially belongs to its own component
        }
    }

    public int getSize() {
        return size;
    }

    public int getEdges() {
        return edges.size();
    }

    public void addEdge(int node1, int node2, int weight) {
        Node getNode1 = nodes.get(node1);
        Node getNode2 = nodes.get(node2);
        Edge edge = new Edge(getNode1, getNode2, weight);
        edges.add(edge);
        getNode1.addNeighbour(edge, getNode2);
        getNode2.addNeighbour(edge, getNode1);

        // Merge components using union-find
        union(node1, node2);
    }

    protected int find(int node) {
        if (parent[node] != node) {
            parent[node] = find(parent[node]); // Path compression
        }
        return parent[node];
    }

    protected void union(int node1, int node2) {
        int root1 = find(node1);
        int root2 = find(node2);
        if (root1 != root2) {
            parent[root1] = root2;
        }
    }

    public boolean isConnected() {
        int root = find(0);
        for (int i = 1; i < size; i++) {
            if (find(i) != root) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String toString() {
        return edges.stream().map(Objects::toString).collect(Collectors.joining(System.lineSeparator()));
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
        for (Node neighbour : current.getNeighbours()) {
            if (!visited.contains(neighbour) && hasCycle(neighbour, target, visited)) {
                return true;
            }
        }
        return false;
    }
}

class KruskalGraph extends AcyclicGraph {
    public KruskalGraph(Graph graph) {
        super(graph.getSize());

        PriorityQueue<Edge> minHeap = new PriorityQueue<>(Comparator.comparingInt(Edge::getWeight));
        minHeap.addAll(graph.edges);

        while (!minHeap.isEmpty() && this.edges.size() < graph.getSize() - 1) {
            Edge minEdge = minHeap.poll();
            int root1 = find(minEdge.getNode(0).getPos());
            int root2 = find(minEdge.getNode(1).getPos());

            if (root1 != root2) {
                super.addEdge(minEdge.getNode(0).getPos(), minEdge.getNode(1).getPos(), minEdge.getWeight());
                union(root1, root2);
            }
        }
    }
}

class PrimGraph extends AcyclicGraph {

    public PrimGraph(Graph graph) {
        super(graph.getSize());
        Set<Node> unvisitedNodes = new HashSet<>(graph.nodes);

        // Use priority queue to select the minimum weight edge efficiently
        PriorityQueue<Edge> minHeap = new PriorityQueue<>(Comparator.comparingInt(Edge::getWeight));
//        Set<Node> visitedNodes = new HashSet<>();

        // Choose a random starting node
        Node startNode = unvisitedNodes.iterator().next();
        unvisitedNodes.remove(startNode);
//        visitedNodes.add(startNode);

        // Add all edges of the start node to the priority queue
        for (Edge edge : startNode.getEdges()) {
            minHeap.offer(edge);
        }

        while (!minHeap.isEmpty()) {

            Edge minEdge = minHeap.poll();
            Node currentNode = minEdge.getNode(0);
            Node neighbourNode = minEdge.getNode(1);

            // If the neighbour node is unvisited, add the edge to the minimum spanning tree
            if (unvisitedNodes.contains(neighbourNode)) {
                long t1 = System.nanoTime();
                super.addEdge(currentNode.getPos(), neighbourNode.getPos(), minEdge.getWeight());
                unvisitedNodes.remove(neighbourNode);
//                visitedNodes.add(neighbourNode);

                // Add all edges of the neighbour node to the priority queue
                for (Edge edge : neighbourNode.getEdges()) {
//                    long t1 = System.nanoTime();
                    if (unvisitedNodes.contains(edge.getNeighbour(neighbourNode))) {
                        minHeap.offer(edge);
                    }
//                    long t2 = System.nanoTime();
//                    time += t2 - t1;
                }
                long t2 = System.nanoTime();
                time += t2 - t1;
            }

        }
    }
}

public class Kruskal {
    public static void main(String[] args) throws IOException {
//        Set<Integer> heap = new HashSet<Integer>();
//        heap.add(98);
//        System.out.println(heap.);
        int upper = 5100;
        Random random = new Random();
        Graph myGraph;
        File file = new File("time4.csv");
        FileWriter fr = new FileWriter(file, false);
        int node1, node2;
        int rep = 3;
        for (int n = 100; n < upper; n += 100) {
            long time, startTime2 = 0, endTime2 = 0;
            long startTime = 0, endTime = 0;
            System.out.println("Graph Size = " + n);
            for (int j = 0; j < rep; j++) {
                myGraph = new Graph(n);
                int i = 0;
//                do {
//                    node1 = random.nextInt(n);
//                    node2 = random.nextInt(n);
//                    if (node1 == node2)
//                        continue;
//                    i++;
//                    myGraph.addEdge(node1, node2, random.nextInt(100));
//                } while (i < n * (n - 1) / 4);
                do {
                    node1 = random.nextInt(n);
                    node2 = random.nextInt(n);
                    if (node1 == node2)
                        continue;
                    myGraph.addEdge(node1, node2, random.nextInt(100));

                } while (!myGraph.isConnected());

                System.out.println("Graph Size = " + myGraph.getEdges());

            startTime += System.nanoTime();
            time = (new PrimGraph(myGraph)).time;
            endTime += System.nanoTime();

//            startTime2 += System.nanoTime();
//            time = (new KruskalGraph(myGraph)).time;
//            endTime2 += System.nanoTime();
        }
            System.out.println("Total-time  = " + (endTime - startTime) / 1e9 / rep);
            fr.write(n + "," + (endTime - startTime) / 1e9 / rep + System.lineSeparator());

//            System.out.println("Total-time  = " + (endTime2 - startTime2) / 1e9 / rep);
//            fr.write(n + "," + (endTime2 - startTime2) / 1e9 / rep + System.lineSeparator());
        }
        fr.close();
    }
}
