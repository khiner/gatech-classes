import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Objects;
import java.util.Comparator;

class Scene {
  Mat4Stack stack;
  float fovDegrees;
  Color backgroundColor;
  List<Light> lights;
  List<Object> objects;

  // Active list of objects to be wrapped in an acceleration data structure.
  // Non-null after `beginAccel` is called, and null again after `endAccel`.
  // When non-null, all objects are added to `accelObjects` instead of `objects`.
  List<Object> accelObjects = null;

  // Named objects are added to the active objects list (`objects` or `accelObjects`)
  // when they are instanced with `createInstance`.
  Map<String, Object> namedObjects;

  Surface surface = null;
  // Active accumulated triangle state:
  Vec3 tri_a = null, tri_b = null, tri_c = null;

  Scene() {
    stack = new Mat4Stack();
    fovDegrees = 0;
    backgroundColor = new Color(0, 0, 0);
    lights = new ArrayList();
    objects = new ArrayList();
    namedObjects = new HashMap();
  }

  void addVertex(Vec3 vertex) {
    vertex = stack.top().transform(vertex);
    if (tri_a == null) tri_a = vertex;
    else if (tri_b == null) tri_b = vertex;
    else if (tri_c == null) tri_c = vertex;
    else throw new IllegalArgumentException("More than three vertices within a single begin/end block.");
  }

  void clearVertices() {
    tri_a = tri_b = tri_c = null;
  }
  void commitVertices() {
    if (surface == null || tri_a == null || tri_b == null || tri_c == null) {
      throw new IllegalArgumentException("Committing vertices without a surface and three points. No triangle added.");
    }

    addObject(new GeometryObject(new Triangle(tri_a, tri_b, tri_c), surface));
  }

  void addBBox(BBox box) {
    final Mat4 transform = stack.top();
    addObject(new GeometryObject(new BBox(transform.transform(box.min), transform.transform(box.max)), surface));
  }

  void nameLatestObject(String name) {
    final Object object = popObject();
    if (object == null) {
      throw new IllegalArgumentException("Attempted to name the latest added object " + name + ", but no " +
        (accelObjects != null ? "accellerated " : "") + "objects have been added to the scene.");
    }

    namedObjects.put(name, object);
  }

  // Create an instance of a named object and add that object to the list of scene objects.
  // Save the inverse of the current transformation matrix as part of the instance.
  void createInstance(String name) {
    if (!namedObjects.containsKey(name)) {
      throw new IllegalArgumentException("Attempted to instance an object named " + name + ", but no object has been given this name.");
    }

    addObject(new InstancedObject(namedObjects.get(name), surface, name, stack.top()));
  }

  void beginAccel() {
    accelObjects = new ArrayList();
  }
  void endAccel() {
    BvhObject bvh = new BvhObject(accelObjects, surface);
    accelObjects = null;
    addObject(bvh);
  }

  // Returns a hit for the intersecting object closest to the ray's origin,
  // or `null` if the ray does not intersect any object.
  Hit raycast(Ray ray) {
    return objects.stream()
      .map(object -> object.raycast(ray))
      .filter(Objects::nonNull)
      .min(Comparator.comparingDouble(hit -> hit.t))
      .orElse(null);
  }

  private void addObject(Object object) {
    (accelObjects != null ? accelObjects : objects).add(object);
  }
  private Object popObject() {
    final List<Object> activeObjects = accelObjects != null ? accelObjects : objects;
    if (activeObjects.isEmpty()) return null;
    return activeObjects.remove(activeObjects.size() - 1);
  }
}
