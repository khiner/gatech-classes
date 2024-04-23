class Sphere extends Geometry {
  final float radius;
  final Vec3 position; // Center

  Sphere(float radius, Vec3 position) {
    this.radius = radius;
    this.position = position;
  }

  Intersection intersect(Ray ray) {
    Vec3 L = position.sub(ray.origin); // Vector from ray origin to sphere center
    float tca = L.dot(ray.direction); // Closest approach along the ray direction
    float d2 = L.dot(L) - tca * tca; // Squared distance from sphere center to the ray at closest approach
    float radius2 = radius * radius;
    if (d2 > radius2) return null; // Ray does not intersect the sphere

    float thc = (float)Math.sqrt(radius2 - d2); // Distance from closest approach to the intersection points
    float t0 = tca - thc;
    float t1 = tca + thc;

    if (t0 > t1) { // Ensure t0 is the smaller value
      float temp = t0;
      t0 = t1;
      t1 = temp;
    }

    if (t0 < 0) {
      t0 = t1; // If t0 is negative, use t1 instead
      if (t0 < 0) return null; // Both t0 and t1 are negative
    }

    Vec3 intersectionPoint = ray.interp(t0);
    Vec3 normal = intersectionPoint.sub(position).normalize();
    return new Intersection(ray, t0, normal);
  }

  BBox getBBox() { return new BBox(position.sub(radius), position.add(radius)); }
}
