// Shadow rays are shot to random locations on this disk.
// When many rays per pixel are used, this creates soft shadows.
class DiskLight extends Light {
  final float radius;
  final Vec3 position, direction;

  DiskLight(float radius, Vec3 position, Vec3 direction, Color c) {
    super(c);
    this.radius = radius;
    this.position = position;
    this.direction = direction.normalize();
  }

  Vec3 samplePoint() {
    // Construct an orthonormal basis `u`/`v` with `direction` as one axis,
    // to rotate the disk to be perpendicular with the `direction`.
    final Vec3 xAxis = new Vec3(1, 0, 0), yAxis = new Vec3(0, 1, 0);
    final Vec3 w = direction; // `w` is aligned with `direction`.
    Vec3 u = xAxis.cross(w); // `u` is perpendicular to both the x-axis and `w`.
    if (u.length() < 1e-6) u = yAxis.cross(w); // If `w` is close to the x-axis, use y-axis instead. 
    u = u.normalize();
    final Vec3 v = w.cross(u);

    final Vec3 p = Random.uniformCircle(radius);
    final Vec3 rotatedP = u.mult(p.x).add(v.mult(p.y)).add(w.mult(p.z));
    return rotatedP.add(position);
  }
}
