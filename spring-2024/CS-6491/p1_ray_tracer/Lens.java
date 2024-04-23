// A lens is a disk that is centered at the origin of the given radius, and that is perpendicular to the z-axis.
// When radius is non-zero, it creates depth of field effects by shooting rays from a lens.
// The origin of the eye rays will not be exactly at (0, 0, 0), but from a random point on the lens.
// `distance` is the distance to the focal plane that is perpendicular to the z-axis.
// All rays for a given pixel, no matter where they originate on the lens, pass through the same point P on the focal plane.
class Lens {
  final float radius, distance;

  Lens(float radius, float distance) {
    this.radius = radius;
    this.distance = distance;
  }

  Vec3 samplePoint() { return radius == 0 ? new Vec3(0, 0, 0) : Random.uniformCircle(radius); }

  // Transform an x/y position on the _image plane_ (defined by `z = imagePoint.z`)
  // into its corresponding point on the _focal plane_.
  Vec3 imageToFocalPoint(Vec3 imagePoint) { return distance == 0 ? imagePoint : imagePoint.mult(distance); }
}
