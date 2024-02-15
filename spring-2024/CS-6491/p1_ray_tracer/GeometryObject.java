class GeometryObject extends Object {
  final Geometry geometry;
  final Surface surface;
  final BBox bbox;

  GeometryObject(Geometry geometry, Surface surface) {
    this.geometry = geometry;
    this.surface = surface;
    this.bbox = geometry.getBBox();
  }

  Hit raycast(Ray ray) {
    final Intersection isect = geometry.intersect(ray);
    if (isect != null) return new Hit(isect, surface);
    return null;
  }

  BBox getBBox() { return bbox; }
}
