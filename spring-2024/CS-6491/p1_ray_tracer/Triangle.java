class Triangle extends Geometry {
  final Vec3 p1, p2, p3;

  Triangle(Vec3 p1, Vec3 p2, Vec3 p3) {
    this.p1 = p1;
    this.p2 = p2;
    this.p3 = p3;
  }

  Vec3 normal() {
    final Vec3 edge1 = p2.sub(p1), edge2 = p3.sub(p1);
    return edge1.cross(edge2).normalize();
  }

  Intersection intersect(Ray ray) {
    final float eps = 1e-5f;

    final Vec3 N = normal();
    // Calculate `t` (the distance from the ray origin to the intersection point).
    final float denom = N.dot(ray.direction);
    if (Math.abs(denom) < eps) return null; // Ray is parallel to the triangle.
  
    final float t = -(N.dot(ray.origin) - N.dot(p1)) / denom;
    if (t < 0) return null; // Intersects behind the ray's origin.

    // Check if the intersection point is inside the triangle.
    final Vec3 p = ray.interp(t);
    final Vec3 edge1 = p2.sub(p1), edge2 = p3.sub(p2), edge3 = p1.sub(p3);
    final boolean
      side1 = N.dot(edge1.cross(p.sub(p1))) > 0,
      side2 = N.dot(edge2.cross(p.sub(p2))) > 0,
      side3 = N.dot(edge3.cross(p.sub(p3))) > 0;
  
    if (side1 && side2 && side3) return new Intersection(ray, t, normal());

    return null;
  }

  BBox getBBox() { return new BBox(Vec3.min(Vec3.min(p1, p2), p3), Vec3.max(Vec3.max(p1, p2), p3)); }
}
