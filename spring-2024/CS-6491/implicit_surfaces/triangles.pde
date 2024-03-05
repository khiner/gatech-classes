// Triangle Mesh

import java.util.List;

List<Vertex> verts;
List<Triangle> triangles;

class Vertex {
  PVector pos;     // position
  PVector normal;  // surface normal
  float r,g,b;     // color

  Vertex (float x, float y, float z) {
    pos = new PVector(x, y, z);
  }
}

class Triangle {
  int v1, v2, v3;
  
  Triangle(int i1, int i2, int i3) {
    v1 = i1;
    v2 = i2;
    v3 = i3;
  }
}

// initialize our list of triangles
void initTriangles() {
  verts = new ArrayList();
  triangles = new ArrayList();
}

// create a new triangle with the given vertex indices
void addTriangle(int i1, int i2, int i3) {
  triangles.add(new Triangle(i1, i2, i3));
}

// add a vertex to the vertex list
int addVertex(PVector p) {
  int index = verts.size();
  Vertex v = new Vertex(p.x, p.y, p.z);
  verts.add (v);
  return (index);
}

// draw the triangles of the surface
void drawSurface() {
  for (int i = 0; i < triangles.size(); i++) {
    Triangle t = triangles.get(i);
    Vertex v1 = verts.get(t.v1), v2 = verts.get(t.v2), v3 = verts.get(t.v3);
    beginShape();
    // add "normal" command before each vertex to use per-vertex (smooth) normals
    vertex(v1.pos.x, v1.pos.y, v1.pos.z);
    vertex(v2.pos.x, v2.pos.y, v2.pos.z);
    vertex(v3.pos.x, v3.pos.y, v3.pos.z);
    endShape(CLOSE);
  }
}

// write triangles to a text file
void writeTriangles(String filename) {
  PrintWriter out = createWriter(filename);
  for (int i = 0; i < triangles.size(); i++) {
    Triangle t = triangles.get(i);
    Vertex v1 = verts.get(t.v1), v2 = verts.get(t.v2), v3 = verts.get(t.v3);
    out.println();
    out.println ("begin");
    out.println ("vertex " + v1.pos.x + " " + v1.pos.y + " " + v1.pos.z);
    out.println ("vertex " + v2.pos.x + " " + v2.pos.y + " " + v2.pos.z);
    out.println ("vertex " + v3.pos.x + " " + v3.pos.y + " " + v3.pos.z);
    out.println ("end");
  }
}
