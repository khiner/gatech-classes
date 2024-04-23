// An `Intersection` represents the intersection of a `Ray` with a `Geometry`.
class Intersection {
  final float t; // Distance along the ray
  final Vec3 point; // Intersection point
  final Vec3 normal; // Normal of the intersected plane

  Intersection(float t, Vec3 point, Vec3 normal) {
    this.t = t;
    this.point = point;
    this.normal = normal;
  }

  Intersection(Ray ray, float t, Vec3 normal) {
    this(t, ray.interp(t), normal.dot(ray.direction) > 0 ? normal.flip() : normal); // Ensure the normal is facing the camera.
  }
}
