class MovingObject extends Object {
  final Object object;
  final Vec3 velocity;

  MovingObject(Object object, Vec3 velocity) {
    this.object = object;
    this.velocity = velocity;
  }

  Hit raycast(Ray ray) {
    final Vec3 translate = velocity.mult(ray.time);
    final Hit hit = object.raycast(Mat4.translate(translate.mult(-1)).transformRay(ray));
    return hit == null ? null : new Hit(hit.t, hit.point.add(translate), hit.normal, hit.surface);
  }

  BBox getBBox() {
    final BBox bbox = object.getBBox();
    return new BBox(bbox.min.sub(velocity), bbox.max.add(velocity)); // Expand bbox by the positive movement
  }
}
