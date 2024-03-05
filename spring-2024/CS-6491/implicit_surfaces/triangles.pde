// Triangle Mesh

import java.util.List;

class Vertex {
  PVector pos;     // position
  PVector normal;  // surface normal
  float r,g,b;     // color

  Vertex(PVector p, PVector n, PVector c) {
    pos = p;
    normal = n;
    r = c.x;
    g = c.y;
    b = c.z;
  }
  color getColor() { return color(r*255, g*255, b*255); }
}

class Triangle {
  int i1, i2, i3;
  
  Triangle(int i1, int i2, int i3) {
    this.i1 = i1;
    this.i2 = i2;
    this.i3 = i3;
  }
}

class Triangles {
  List<Vertex> verts;
  List<Triangle> triangles;

  Triangles() {
    verts = new ArrayList();
    triangles = new ArrayList();
  }

  int triangleCount() { return triangles.size(); }
  int vertexCount() { return verts.size(); }

  void add(int i1, int i2, int i3) { triangles.add(new Triangle(i1, i2, i3)); }

  int addVertex(PVector p, PVector n, PVector c) {
    verts.add(new Vertex(p, n, c));
    return verts.size() - 1;
  }

  void draw() {
    for (Triangle t : triangles) {
      final Vertex v1 = verts.get(t.i1), v2 = verts.get(t.i2), v3 = verts.get(t.i3);
      beginShape();
      if (draw_flags.smooth_normals) normal(v1.normal.x, v1.normal.y, v1.normal.z);
      fill(v1.getColor());
      vertex(v1.pos.x, v1.pos.y, v1.pos.z);
      if (draw_flags.smooth_normals) normal(v2.normal.x, v2.normal.y, v2.normal.z);
      fill(v2.getColor());
      vertex(v2.pos.x, v2.pos.y, v2.pos.z);
      if (draw_flags.smooth_normals) normal(v3.normal.x, v3.normal.y, v3.normal.z);
      fill(v3.getColor());
      vertex(v3.pos.x, v3.pos.y, v3.pos.z);
      endShape(CLOSE);
    }
  }

  // write triangles to a text file
  void write(String filename) {
    PrintWriter out = createWriter(filename);
    for (Triangle t : triangles) {
      final Vertex v1 = verts.get(t.i1), v2 = verts.get(t.i2), v3 = verts.get(t.i3);
      out.println();
      out.println("begin");
      out.println("vertex " + v1.pos.x + " " + v1.pos.y + " " + v1.pos.z);
      out.println("vertex " + v2.pos.x + " " + v2.pos.y + " " + v2.pos.z);
      out.println("vertex " + v3.pos.x + " " + v3.pos.y + " " + v3.pos.z);
      out.println("end");
    }
    out.close();
  }
}
