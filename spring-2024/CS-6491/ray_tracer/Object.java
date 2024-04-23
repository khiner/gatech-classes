// All objects rendered in a `Scene` inherit from `Object`.
abstract class Object {
  abstract Hit raycast(Ray ray);
  abstract BBox getBBox();
}
