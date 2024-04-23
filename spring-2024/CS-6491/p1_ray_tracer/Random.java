class Random {
  // Generate a point on the unit disk, with z = 0.
  // It's easy enough to generate a proper uniform distribution within a circle,
  // so I'm using that rather than rejection sampling for efficiency.
  static Vec3 uniformCircle(float radius) {
    final float theta = (float)(2 * Math.PI * Math.random());
    final float r = (float)(radius * Math.sqrt(Math.random()));
    return new Vec3((float)Math.cos(theta), (float)Math.sin(theta), 0).mult(r);
  }

  // Generate a random 3D point inside a sphere of a given radius, using rejection sampling.
  static Vec3 uniformSphere(float radius) {
    Vec3 p;
    do {
      p = new Vec3((float)Math.random(), (float)Math.random(), (float)Math.random()).mult(2).sub(1);
    } while (p.dot(p) >= 1.0); // Ensure the point is within the unit sphere.

    return p.mult(radius);
  }
}
