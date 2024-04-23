class InstancedObject extends Object {
  final Surface surface; // Overrides the object's surface if not null
  final Object object;
  final String name;
  final Mat4 transform, invTransform;

  InstancedObject(Object object, Surface surface, String name, Mat4 transform) {
    this.object = object;
    this.surface = surface;
    this.name = name;
    this.transform = transform;
    this.invTransform = transform.invert();
    if (invTransform == null) throw new RuntimeException("Attempted to instance an object with an uninvertible transform.");
  }

  Hit raycast(Ray ray) {
    final Hit hit = object.raycast(invTransform.transformRay(ray));
    if (hit != null) {
      return new Hit(hit.t,
        transform.transform(hit.point), // Transform the point back to world space.
        invTransform.transpose().transformDirection(hit.normal).normalize(),
        surface != null ? surface : hit.surface
        );
    }
    return hit;
  }

  BBox getBBox() {
    final BBox bbox = object.getBBox();
    return new BBox(transform.transform(bbox.min), transform.transform(bbox.max));
  }
}
