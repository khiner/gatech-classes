// A `Hit` is an `Intersection` with rendering properties (currently just a surface).
class Hit extends Intersection {
  final Surface surface;

  Hit(float t, Vec3 point, Vec3 normal, Surface surface) {
    super(t, point, normal);
    this.surface = surface;
  }
  Hit(Ray ray, float t, Vec3 normal, Surface surface) {
    super(ray, t, normal);
    this.surface = surface;
  }
  Hit(Intersection isect, Surface surface) {
    this(isect.t, isect.point, isect.normal, surface);
  }
}
