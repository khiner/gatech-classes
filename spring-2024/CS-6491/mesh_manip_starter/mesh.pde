import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Random;
import java.util.Objects;
import java.util.function.BiConsumer;

import processing.core.PVector;

// Avoid floating-point precision issues when using `PVector`s as keys.
String toKey(PVector v) { return String.format("%1.5f,%1.5f,%1.5f", v.x, v.y, v.z); }

Random rand = new Random();

class Vertex {
  PVector position;
  PVector normal;
  HalfEdge edge; // One outgoing half-edge

  Vertex(PVector position) { this.position = position; }
  Vertex(float x, float y, float z) { this(new PVector(x, y, z)); }
}

class Face {
  HalfEdge edge; // One half-edge in the face
  PVector normal;
  color col = color(255, 255, 255);
  int numVertices = 0;

  PVector calculateCentroid() {
    PVector centroid = new PVector();
    HalfEdge e = edge;
    do { centroid.add(e.target.position); } while ((e = e.next) != edge);
    return centroid.div(numVertices);
  }
}

class HalfEdge {
  final Vertex target; // The vertex the edge points to
  HalfEdge next, opposite;
  Face face; // The face left of the edge

  HalfEdge(Vertex target) { this.target = target; }
}

// Convenience helpers
PVector vec(float x, float y) { return new PVector(x, y); }
PVector vec(float x, float y, float z) { return new PVector(x, y, z); }

// Based on https://stackoverflow.com/a/7898685/780425
PVector hsvToRgb(float hue, float saturation, float v) {
    int h = (int)(hue*6);
    float f = hue*6 - h;
    float p = v*(1 - saturation);
    float q = v*(1 - f*saturation);
    float t = v*(1 - (1 - f)*saturation);
    switch (h) {
      case 0: return vec(v, t, p);
      case 1: return vec(q, v, p);
      case 2: return vec(p, v, t);
      case 3: return vec(p, q, v);
      case 4: return vec(t, p, v);
      case 5: return vec(v, p, q);
      default: throw new RuntimeException("Something went wrong when converting from HSV to RGB. Input was " + hue + ", " + saturation + ", " + v);
    }
}

color randCol() {
  PVector v = hsvToRgb(rand.nextFloat(), 0.4, 1.0);
  return color(v.x * 255, v.y * 255, v.z * 255);
}

enum EdgeMove { Next, Previous, Opposite, Swing }

class Mesh {
  List<Vertex> vertices = new ArrayList<>();
  List<Face> faces = new ArrayList<>();
  List<HalfEdge> edges = new ArrayList<>();

  Mesh(List<PVector> vertices, List<int[]> faces) {
    for (PVector v : vertices) this.vertices.add(new Vertex(v));
    for (int[] face : faces) addFace(face);
    setOpposites();
    calculateNormals();
  }

  Mesh(List<PVector> vertices, List<int[]> faces, boolean smoothShading, boolean randomColors, boolean renderEdges) {
    this(vertices, faces);
    this.smoothShading = smoothShading;
    this.randomColors = randomColors;
    this.renderEdges = renderEdges;
    updateColors();
  }

  Mesh(Mesh other, List<PVector> vertices, List<int[]> faces) {
    this(vertices, faces, other.smoothShading, other.randomColors, other.renderEdges);
  }

  boolean smoothShading = false, randomColors = false, renderEdges = false;

  HalfEdge debugEdge = null; // Current visualized edge.
  boolean showDebugEdge = false; // Visualize the `currEdge`

  void updateColors() {
    for (Face face : faces) face.col = randomColors ? randCol() : color(255, 255, 255);
  }

  void toggleSmoothShading() { smoothShading = !smoothShading; }
  void toggleRandomColors() {
    randomColors = !randomColors;
    updateColors();
  }
  void toggleEdges() { renderEdges = !renderEdges; }
  void toggleEdgeDebug() {
    showDebugEdge = !showDebugEdge;
    if (showDebugEdge && debugEdge == null && !edges.isEmpty()) debugEdge = edges.get(0);
  }
  void moveDebugEdge(EdgeMove move) { debugEdge = getEdge(debugEdge, move); }

  void addFace(int[] vertexIndices) {
    Face face = new Face();
    HalfEdge firstEdge = null, prevEdge = null;
    for (int i = 0; i < vertexIndices.length; i++) {
      Vertex currentVertex = vertices.get(vertexIndices[i]);
      HalfEdge edge = new HalfEdge(currentVertex);
      if (i == 0) firstEdge = edge;
      else prevEdge.next = edge;

      edge.face = face;
      if (currentVertex.edge == null) currentVertex.edge = edge;

      edges.add(edge);
      prevEdge = edge;
    }

    prevEdge.next = firstEdge; // Close the loop.
    face.edge = firstEdge;
    face.numVertices = vertexIndices.length;
    faces.add(face);
  }

  void setOpposites() {
    class EdgeKey {
      PVector start, end;

      EdgeKey(PVector start, PVector end) {
        this.start = start;
        this.end = end;
      }

      @Override
      public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof EdgeKey)) return false;
        EdgeKey edgeKey = (EdgeKey)o;
        return start.equals(edgeKey.start) && end.equals(edgeKey.end);
      }
      @Override
      public int hashCode() {
        return Objects.hash(start.x, start.y, start.z, end.x, end.y, end.z);
      }
    }

    Map<EdgeKey, HalfEdge> edgeMap = new HashMap<>();
    for (HalfEdge edge : edges) {
      EdgeKey edgeKey = new EdgeKey(edge.target.position, edge.next.target.position);
      EdgeKey oppositeKey = new EdgeKey(edge.next.target.position, edge.target.position);
      if (edgeMap.containsKey(oppositeKey)) {
        HalfEdge oppositeEdge = edgeMap.get(oppositeKey);
        edge.opposite = oppositeEdge;
        oppositeEdge.opposite = edge;
      } else {
        edgeMap.put(edgeKey, edge);
      }
    }
  }

  void calculateNormals() {
    for (Vertex vertex : vertices) vertex.normal = new PVector();

    for (Face face : faces) {
      HalfEdge edge1 = face.edge, edge2 = edge1.next, edge3 = edge2.next;
      final PVector v1 = PVector.sub(edge2.target.position, edge1.target.position);
      final PVector v2 = PVector.sub(edge3.target.position, edge1.target.position);
      face.normal = v1.cross(v2).normalize();

      // Accumulate vertex normals (we normalize later).
      do {
        edge1.target.normal.add(face.normal);
        edge1 = edge1.next;
      } while (edge1 != face.edge);
    }

    for (Vertex vertex : vertices) vertex.normal.normalize();
  }

  void addRandomNoise() {
    for (Vertex vertex : vertices) {
      // Add random value to each vertex in range [-.1, .1] * surface normal.
      vertex.position.add(PVector.mult(vertex.normal, (2*rand.nextFloat() - 1) * 0.1));
    }
    calculateNormals();
  }

  void smoothLaplacian(float lambda) {
    List<PVector> newPositions = new ArrayList<>(vertices.size());
    for (Vertex vertex : vertices) newPositions.add(vertex.position.copy());

    for (int i = 0; i < vertices.size(); i++) {
      Vertex v = vertices.get(i);
      List<Vertex> neighbors = getVertexNeighbors(v);

      PVector centroid = new PVector(0, 0, 0);
      for (Vertex neighbor : neighbors) centroid.add(neighbor.position);
      if (!neighbors.isEmpty()) centroid.div(neighbors.size());

      PVector originalPosition = v.position;
      PVector movementVector = PVector.sub(centroid, originalPosition).mult(lambda);
      newPositions.set(i, PVector.add(originalPosition, movementVector));
    }

    for (int i = 0; i < vertices.size(); i++) vertices.get(i).position = newPositions.get(i);

    calculateNormals();
  }

  void smoothTaubin(float lambda1, float mu) {
    smoothLaplacian(lambda1);
    smoothLaplacian(mu);
  }

  // Find the neighboring vertices of a vertex.
  List<Vertex> getVertexNeighbors(Vertex vertex) {
    List<Vertex> neighbors = new ArrayList<>();
    HalfEdge startEdge = vertex.edge, edge = startEdge;
    do {
      if (edge.opposite != null) neighbors.add(edge.opposite.target);
    } while ((edge = getSwingEdge(edge)) != null && edge != startEdge);
    return neighbors;
  }
  List<Face> getAdjacentFaces(Vertex vertex) {
    List<Face> adjacentFaces = new ArrayList<>();
    HalfEdge startEdge = vertex.edge, edge = startEdge;
    do { adjacentFaces.add(edge.face); } while ((edge = getSwingEdge(edge)) != startEdge);
    return adjacentFaces;
  }

  HalfEdge getEdge(HalfEdge edge, EdgeMove move) {
    if (edge == null) return null;

    switch (move) {
      case Next: return getNextEdge(edge);
      case Previous: return getPreviousEdge(edge);
      case Opposite: return getOppositeEdge(edge);
      case Swing: return getSwingEdge(edge);
      default: return null;
    }
  }

  HalfEdge getNextEdge(HalfEdge edge) { return edge.next; }
  HalfEdge getPreviousEdge(HalfEdge currentEdge) {
    // Iterate over face edges to find the one whose `next is `currentEdge`.
    HalfEdge edge = currentEdge;
    while (edge.next != currentEdge) edge = edge.next;
    return edge;
  }
  HalfEdge getOppositeEdge(HalfEdge edge) { return edge.opposite; }
  HalfEdge getSwingEdge(HalfEdge edge) { return edge.opposite != null ? edge.opposite.next : null; }

  void visualizeDirectedEdge(HalfEdge edge) {
    if (edge == null || !showDebugEdge) return;

    final int numSpheres = 3;
    final int sphereSteps = 15;
    final float sphereStepRatio = 1.f/sphereSteps;

    final Vertex start = edge.target, end = edge.next.target;
    final PVector
      edgeVec = PVector.sub(end.position, start.position),
      sphereStep = PVector.mult(edgeVec, sphereStepRatio),
      sphereStart = PVector.add(start.position, PVector.mult(sphereStep, (1 + sphereSteps - numSpheres)/2.f));
    final float
      sphereStepMag = sphereStep.mag(),
      sphereRadMin = sphereStepMag/2,
      sphereRadMax = sphereStepMag, sphereRadInc = (sphereRadMax - sphereRadMin)/(numSpheres - 1);
    final PVector
      sphereOffset = PVector.cross(edge.face.normal, edgeVec.copy().normalize(), null).normalize().mult(sphereStepMag);

    fill(255, 100, 100);
    // Draw spheres along the edge.
    for (int i = 0; i < numSpheres; i++) {
      PVector pos = PVector.add(sphereStart, PVector.mult(sphereStep, i)).add(sphereOffset);
      pushMatrix();
      translate(pos.x, pos.y, pos.z);
      sphere(sphereRadMin + (numSpheres - i - 1)*sphereRadInc);
      popMatrix();
    }
  }

  Mesh createDual() {
    Map<Face, Integer> faceToVertexIndex = new HashMap<>(); // Original face -> corresponding dual vertex

    List<PVector> dualVertices = new ArrayList<>(vertices.size());
    List<int[]> dualFaces = new ArrayList<>(faces.size());

    // Create a vertex for each face in the original mesh.
    for (Face face : faces) {
      faceToVertexIndex.put(face, dualVertices.size());
      dualVertices.add(face.calculateCentroid());
    }
    // Create a face for each vertex in the original mesh.
    for (Vertex vertex : vertices) {
      List<Face> adjacentFaces = getAdjacentFaces(vertex);
      int[] dualFaceVertices = new int[adjacentFaces.size()];
      for (int i = 0; i < adjacentFaces.size(); i++) dualFaceVertices[i] = faceToVertexIndex.get(adjacentFaces.get(i));

      dualFaces.add(dualFaceVertices);
    }

    return new Mesh(this, dualVertices, dualFaces);
  }

  Mesh subdivideMidpoint() {
    List<PVector> newVertices = new ArrayList<>();
    Map<String, Integer> indexForVertex = new HashMap<>();
    List<int[]> newFaces = new ArrayList<>();    
    for (Face f : faces) {
      int[] midpointsIndices = new int[f.numVertices];
      HalfEdge edge = f.edge;
      int prevMidIndex = -1, firstMidIndex = -1;
      for (int i = 0; i < f.numVertices; i++) {
        // `normalize` to project onto unit sphere
        final PVector startPos = edge.target.position.normalize();
        final PVector endPos = edge.next.target.position.normalize();
        final PVector midPos = PVector.add(startPos, endPos).div(2).normalize();
        if (indexForVertex.putIfAbsent(toKey(startPos), newVertices.size()) == null) newVertices.add(startPos);
        if (indexForVertex.putIfAbsent(toKey(midPos), newVertices.size()) == null) newVertices.add(midPos);

        final int midIndex = indexForVertex.get(toKey(midPos));
        midpointsIndices[i] = midIndex;

        if (prevMidIndex != -1) newFaces.add(new int[]{indexForVertex.get(toKey(startPos)), midIndex, prevMidIndex});
        else firstMidIndex = midIndex; // Save the first midpoint to close the loop later.

        prevMidIndex = midIndex;
        edge = edge.next;
      }

      // Close the loop, and add face containing all midpoints.
      newFaces.add(new int[]{indexForVertex.get(toKey(f.edge.target.position.normalize())), firstMidIndex, prevMidIndex});
      newFaces.add(midpointsIndices);
    }

    return new Mesh(this, newVertices, newFaces);
  }

  Mesh subdivideCatmullClark() {
    Map<Face, PVector> facePoints = new HashMap<>();
    for (Face face : this.faces) facePoints.put(face, face.calculateCentroid());

    Map<HalfEdge, PVector> edgePoints = new HashMap<>();
    for (Face face : this.faces) {
      HalfEdge edge = face.edge;
      do {
        if (!edgePoints.containsKey(edge)) { // Avoid duplicating calculations for shared edges
          PVector edgePoint = PVector.add(edge.target.position, edge.next.target.position)
                                .add(facePoints.get(face))
                                .add(facePoints.get(edge.opposite.face))
                                .div(4);
          edgePoints.put(edge, edgePoint);
          edgePoints.put(edge.opposite, edgePoint); // Shared edge, same point
        }
      } while ((edge = edge.next) != face.edge);
    }

    Map<Vertex, PVector> newPositions = new HashMap<>();
    for (Vertex vertex : this.vertices) {
      PVector avgEdgePoint = new PVector(), avgFacePoint = new PVector();
      HalfEdge edge = vertex.edge;
      int count = 0;
      do {
        avgEdgePoint.add(edgePoints.get(edge));
        avgFacePoint.add(facePoints.get(edge.face));
        count++;
      } while ((edge = getSwingEdge(edge)) != vertex.edge);

      avgEdgePoint.div(count).mult(2);
      avgFacePoint.div(count);
      newPositions.put(vertex, PVector.add(avgEdgePoint, avgFacePoint).add(vertex.position.copy().mult(count - 3)).div(count));
    }

    // Construct the final mesh.
    List<PVector> newVertices = new ArrayList<>();

    // Cached vertex add.
    // Note: Unlike in midpoint subdivision, we don't seem to use a rounded string key here.
    Map<PVector, Integer> indexForVertex = new HashMap<>();
    BiConsumer<PVector, List<Integer>> addVertex = (vertex, currFace) -> {
      if (indexForVertex.putIfAbsent(vertex, newVertices.size()) == null) newVertices.add(vertex);
      currFace.add(indexForVertex.get(vertex));
    };

    List<int[]> newFaces = new ArrayList<>();
    for (Face face : this.faces) {
      HalfEdge edge = face.edge;
      do {
        List<Integer> fis = new ArrayList<>(4); // Face indices
        addVertex.accept(newPositions.get(edge.target), fis);
        addVertex.accept(edgePoints.get(edge), fis);
        addVertex.accept(facePoints.get(face), fis);
        addVertex.accept(edgePoints.get(getPreviousEdge(edge)), fis);
        newFaces.add(fis.stream().mapToInt(Integer::intValue).toArray());
      } while ((edge = edge.next) != face.edge);
    }

    return new Mesh(this, newVertices, newFaces);
  }

  void draw() {
    for (Face face : faces) {
      fill(face.col);
      if (renderEdges) stroke(0); // Black edges
      else noStroke();

      beginShape();
      HalfEdge edge = face.edge;
      do {
        if (smoothShading) normal(edge.target.normal.x, edge.target.normal.y, edge.target.normal.z);
        else normal(face.normal.x, face.normal.y, face.normal.z);

        Vertex v = edge.target;
        vertex(v.position.x, v.position.y, v.position.z);
        edge = edge.next;
      } while (edge != face.edge);
      endShape(CLOSE);
    }
    
    if (debugEdge != null && showDebugEdge) {
      visualizeDirectedEdge(debugEdge);
    }
  }
}

Mesh loadMesh(String filename) {
  final String[] lines = loadStrings(filename);
  final int numVertices = int(split(lines[0], " ")[1]);
  final int numFaces = int(split(lines[1], " ")[1]);

  List<PVector> vertices = new ArrayList<>(numVertices);
  List<int[]> faces = new ArrayList<>(numFaces);

  // Read vertices.
  String[] words;
  for (int i = 0; i < numVertices; i++) {
    words = split(lines[i + 2], " ");
    vertices.add(new PVector(float(words[0]), float(words[1]), float(words[2])));
  }
  
  // Read faces.
  for (int i = 0; i < numFaces; i++) {
    words = split(lines[i + numVertices + 2], " ");
    final int nVerts = int(words[0]);
    final int[] vertexIndices = new int[nVerts];
    for (int k = 1; k <= nVerts; k++) vertexIndices[k - 1] = int(words[k]);
    faces.add(vertexIndices);
  }

  Mesh mesh = new Mesh(vertices, faces);
  println("Mesh loaded: " + mesh.vertices.size() + " vertices, " + mesh.faces.size() + " faces.");
  return mesh;
}
