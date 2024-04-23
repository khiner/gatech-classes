import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Objects;
import java.util.Comparator;

class Scene {
  final Vec3 cameraPosition = new Vec3(0, 0, 0); // Not yet configurable
  Lens lens = new Lens(0, 0);
  float fovDegrees = 0;
  // When the number of rays is larger than one, sub-pixel rays are created in random positions within each pixel.
  // The colors of these rays are averaged together to give the final color of the pixel.
  int raysPerPixel = 1;
  Color backgroundColor = new Color(0, 0, 0);
  List<Light> lights = new ArrayList();
  List<Object> objects = new ArrayList();

  Mat4Stack stack = new Mat4Stack();

  // Active list of objects to be wrapped in an acceleration data structure.
  // Non-null after `beginAccel` is called, and null again after `endAccel`.
  // When `accelObjects != null`, all objects are added to `accelObjects` instead of `objects`.
  List<Object> accelObjects = null;

  // Named objects are added to the active objects list (`objects` or `accelObjects`)
  // when they are instanced with `createInstance`.
  Map<String, Object> namedObjects = new HashMap();

  Surface surface = null;
  // Active accumulated triangle state:
  int tri_i = 0; // Current triangle vertex index (0/1/2)
  Vec3[] tri_verts = {null, null, null};

  void addVertex(Vec3 vertex) {
    tri_i = (tri_i + 1) % 3;
    tri_verts[tri_i] = stack.top().transform(vertex);
  }

  void clearVertices() { tri_i = -1; }
  void commitVertices() {
    if (tri_i == -1 || tri_i > 2) throw new IllegalArgumentException("Committing vertices without adding three vertices first. No triangle added.");

    addObject(new GeometryObject(new Triangle(tri_verts[0], tri_verts[1], tri_verts[2]), surface));
  }

  void addBBox(BBox box) {
    final Mat4 transform = stack.top();
    addObject(new GeometryObject(new BBox(transform.transform(box.min), transform.transform(box.max)), surface));
  }
  
  void addSphere(Sphere sphere) {
    final Mat4 transform = stack.top();
    addObject(new GeometryObject(new Sphere(sphere.radius, transform.transform(sphere.position)), surface));
  }

  void nameLatestObject(String name) {
    final Object object = popObject();
    if (object == null) {
      throw new IllegalArgumentException("Attempted to name the latest added object " + name + ", but no " +
        (accelObjects != null ? "accellerated " : "") + "objects have been added to the scene.");
    }

    namedObjects.put(name, object);
  }

  // `velocity` is the amount of translation the object undergoes during one frame.
  void setLatestObjectVelocity(Vec3 velocity) {
    final Object object = popObject();
    if (object == null) {
      throw new IllegalArgumentException("Attempted to set the velocity of the latest added object, but no " +
        (accelObjects != null ? "accellerated " : "") + "objects have been added to the scene.");
    }
    addObject(new MovingObject(object, velocity));
  }

  // Create an instance of a named object and add that object to the list of scene objects.
  // Save the inverse of the current transformation matrix as part of the instance.
  void createInstance(String name) {
    if (!namedObjects.containsKey(name)) {
      throw new IllegalArgumentException("Attempted to instance an object named " + name + ", but no object has been given this name.");
    }

    addObject(new InstancedObject(namedObjects.get(name), surface, name, stack.top()));
  }

  void beginAccel() { accelObjects = new ArrayList(); }
  void endAccel() {
    BvhObject bvh = new BvhObject(accelObjects, surface);
    accelObjects = null;
    addObject(bvh);
  }

  // Returns a hit for the intersecting object closest to the ray's origin,
  // or `null` if the ray does not intersect any object.
  Hit raycast(Ray ray) {
    return objects.stream()
      .map(o -> o.raycast(ray))
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
