class DiskLight extends Light {
  final float radius;
  final Vec3 position, direction;

  DiskLight(float radius, Vec3 position, Vec3 direction, Color c) {
    super(c);
    this.radius = radius;
    this.position = position;
    this.direction = direction.normalize();
  }

  Vec3 samplePosition() {
    // Generate a point on the unit disk
    final float theta = (float)(2 * Math.PI * Math.random());
    final float r = (float)(radius * Math.sqrt(Math.random()));
    final Vec3 diskPoint = new Vec3((float)Math.cos(theta), (float)Math.sin(theta), 0).mult(r);

    // Construct an orthonormal basis `u`/`v` with 'direction' as one axis.
    final Vec3 xAxis = new Vec3(1, 0, 0), yAxis = new Vec3(0, 1, 0); // Standard basis vectors
    final Vec3 w = direction; // 'w' is aligned with the 'direction' of the disk.
    Vec3 u = xAxis.cross(w); // 'u' is perpendicular to both the x-axis and 'w'.
    if (u.length() < 1e-6) u = yAxis.cross(w); // Use y-axis instead if 'w' is close to the x-axis.
    u = u.normalize();
    final Vec3 v = w.cross(u);

    final Vec3 rotatedPoint = u.mult(diskPoint.x).add(v.mult(diskPoint.y)).add(w.mult(diskPoint.z));
    return rotatedPoint.add(position);
  }
}
