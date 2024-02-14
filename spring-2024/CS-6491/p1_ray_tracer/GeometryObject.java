class GeometryObject extends Object {
  final Geometry geometry;
  final Surface surface;

  GeometryObject(Geometry geometry, Surface surface) {
    this.geometry = geometry;
    this.surface = surface;
  }

  Hit raycast(Ray ray) {
    final Intersection isect = geometry.intersect(ray);
    if (isect != null) return new Hit(isect, surface);
    return null;
  }

  BBox getBBox() {
    return geometry.getBBox();
  }
}
