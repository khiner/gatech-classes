import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Random;

import processing.core.PVector;

Random rand = new Random();

class Vertex {
  final PVector position;
  PVector normal;
  HalfEdge halfEdge; // One outgoing half-edge

  Vertex(float x, float y, float z) {
    this.position = new PVector(x, y, z);
  }
}

class Face {
  HalfEdge halfEdge; // One half-edge in the face
  PVector normal;
  color col = color(255, 255, 255);
}

class HalfEdge {
  final Vertex target; // The vertex the edge points to
  HalfEdge next;
  HalfEdge opposite;
  Face face; // The face left of the edge

  HalfEdge(Vertex target) {
    this.target = target;
  }
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


class Mesh {
  List<Vertex> vertices = new ArrayList<>();
  List<Face> faces = new ArrayList<>();
  List<HalfEdge> halfEdges = new ArrayList<>();

  boolean smoothShading = false, randomColors = false, renderEdges = false;

  void toggleSmoothShading() { smoothShading = !smoothShading; }
  void toggleRandomColors() {
    randomColors = !randomColors;
    for (Face face : faces) {
      face.col = randomColors ? randCol()  : color(255, 255, 255);
    }
  }
  void toggleEdges() { renderEdges = !renderEdges; }

  void addFace(int[] vertexIndices) {
    Face face = new Face();
    HalfEdge firstHalfEdge = null, prevHalfEdge = null;
    for (int i = 0; i < vertexIndices.length; i++) {
      Vertex currentVertex = vertices.get(vertexIndices[i]);
      HalfEdge halfEdge = new HalfEdge(currentVertex);
      if (i == 0) firstHalfEdge = halfEdge;
      else prevHalfEdge.next = halfEdge;
      
      halfEdge.face = face;
      if (currentVertex.halfEdge == null) currentVertex.halfEdge = halfEdge;
      
      halfEdges.add(halfEdge);
      prevHalfEdge = halfEdge;
    }

    prevHalfEdge.next = firstHalfEdge; // Close the loop.
    face.halfEdge = firstHalfEdge;
    faces.add(face);
  }
  
  // Set opposite half-edges
  void setOpposites() {
    Map<String, HalfEdge> edgeMap = new HashMap<>();
    for (HalfEdge he : halfEdges) {
      final String edgeKey = he.target.position + " " + he.next.target.position;
      final String oppositeKey = he.next.target.position + " " + he.target.position;
      if (edgeMap.containsKey(oppositeKey)) {
        HalfEdge oppositeHalfEdge = edgeMap.get(oppositeKey);
        he.opposite = oppositeHalfEdge;
        oppositeHalfEdge.opposite = he;
      } else {
        edgeMap.put(edgeKey, he);
      }
    }
  }

  void calculateNormals() {
    for (Vertex vertex : vertices) vertex.normal = new PVector();

    for (Face face : faces) {
      HalfEdge edge1 = face.halfEdge, edge2 = edge1.next, edge3 = edge2.next;
      final PVector v1 = PVector.sub(edge1.target.position, edge2.target.position);
      final PVector v2 = PVector.sub(edge1.target.position, edge3.target.position);
      face.normal = v2.cross(v1).normalize();

      // Accumulate vertex normals (we normalize later).
      do {
        edge1.target.normal.add(face.normal);
        edge1 = edge1.next;
      } while (edge1 != face.halfEdge);
    }

    // Normalize vertex normals
    for (Vertex vertex : vertices) vertex.normal.normalize();
  }

  void draw() {
    for (Face face : faces) {
      fill(face.col);
      if (renderEdges) stroke(0); // Black edges
      else noStroke();

      beginShape();
      HalfEdge edge = face.halfEdge;
      do {
        if (smoothShading) {
          normal(edge.target.normal.x, edge.target.normal.y, edge.target.normal.z);
        } else {
          normal(face.normal.x, face.normal.y, face.normal.z);
        }

        Vertex v = edge.target;
        vertex(v.position.x, v.position.y, v.position.z);
        edge = edge.next;
      } while (edge != face.halfEdge);
      endShape(CLOSE);
    }
  }
}

Mesh loadMesh(String filename) {
  Mesh mesh = new Mesh();
  final String[] lines = loadStrings(filename);
  final int numVertices = int(split(lines[0], " ")[1]);
  final int numFaces = int(split(lines[1], " ")[1]);

  // Read vertices.
  String[] words;
  for (int i = 0; i < numVertices; i++) {
    words = split(lines[i + 2], " ");
    mesh.vertices.add(new Vertex(float(words[0]), float(words[1]), float(words[2])));
  }
  
  // Read faces.
  for (int i = 0; i < numFaces; i++) {
    words = split(lines[i + numVertices + 2], " ");

    final int nVerts = int(words[0]);
    int[] vertexIndices = new int[nVerts];
    for (int k = 1; k <= nVerts; k++) vertexIndices[k - 1] = int(words[k]);
    
    mesh.addFace(vertexIndices);
  }
  
  mesh.setOpposites(); // Set opposite half-edges.
  mesh.calculateNormals();
  
  println("Mesh loaded: " + mesh.vertices.size() + " vertices, " + mesh.faces.size() + " faces.");
  return mesh;
}
