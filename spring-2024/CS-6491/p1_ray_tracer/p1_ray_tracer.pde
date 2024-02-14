// This is the starter code for the CS 6491 Ray Tracing project.
//
// The most important part of the code is the interpreter, which will help
// you parse the scene description (.cli) files.

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

void setup() {
  size (300, 300);
  noStroke();
  background (0, 0, 0);
}

static Integer charToCliFileNumber(char ch) {
  if (ch >= '1' && ch <= '9') return ch - '0';
  if (ch == '0') return 10;
  if (ch == 'a') return 11;
  return null;
}

void keyPressed() {
  final Integer cli_file_number = charToCliFileNumber(key);
  if (cli_file_number == null) return;

  interpret(String.format("s%02d.cli", cli_file_number), new Scene());
}

/**** Scene description ****/

class Color {
  final float r, g, b; // In range [0, 1]

  Color(float r, float g, float b) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  
  Color add(Color other) { return new Color(r + other.r, g + other.g, b + other.b); }
  Color add(float scalar) { return new Color(r + scalar, g + scalar, b + scalar); }

  Color mult(Color other) { return new Color(r * other.r, g * other.g, b * other.b); }
  Color mult(float scalar) { return new Color(r * scalar, g * scalar, b * scalar); }

  color get() { return color(r * 255, g * 255, b * 255); }
}

class Light {
  final Vec3 position;
  final Color c;

  Light(Vec3 position, Color c) {
    this.position = position;
    this.c = c;
  }
}

class Surface {
  final Color diffuse;

  Surface(Color diffuse) {
    this.diffuse = diffuse;
  }
}

// An `Intersection` represents the intersection of a `Ray` with a `Geometry`.
class Intersection {
  final float t; // Distance along the ray
  final Vec3 point; // Intersection point
  final Vec3 normal; // Normal of the intersected plane

  Intersection(float t, Vec3 point, Vec3 normal) {
    this.t = t;
    this.point = point;
    this.normal = normal;
  }
  Intersection(Ray ray, float t, Vec3 normal) {
    this(t, ray.interp(t), normal.dot(ray.direction) > 0 ? normal.flip() : normal); // Ensure the normal is facing the camera.
  }
}

abstract class Geometry {
  abstract Intersection intersect(Ray ray);
}

class Triangle extends Geometry {
  final Vec3 p1, p2, p3;

  Triangle(Vec3 p1, Vec3 p2, Vec3 p3) {
    this.p1 = p1;
    this.p2 = p2;
    this.p3 = p3;
  }

  Vec3 normal() {
    final Vec3 edge1 = p2.sub(p1), edge2 = p3.sub(p1);
    return edge1.cross(edge2).normalize();
  }

  Intersection intersect(Ray ray) {
    final float eps = 1e-5;

    final Vec3 N = normal();
    // Calculate `t` (the distance from the ray origin to the intersection point).
    final float denom = N.dot(ray.direction);
    if (abs(denom) < eps) return null; // Ray is parallel to the triangle.
  
    final float t = -(N.dot(ray.origin) - N.dot(p1)) / denom;
    if (t < 0) return null; // Intersects behind the ray's origin.

    // Check if the intersection point is inside the triangle.
    final Vec3 p = ray.interp(t);
    final Vec3 edge1 = p2.sub(p1), edge2 = p3.sub(p2), edge3 = p1.sub(p3);
    final boolean
      side1 = -N.dot(edge1.cross(p.sub(p1))) < eps,
      side2 = -N.dot(edge2.cross(p.sub(p2))) < eps,
      side3 = -N.dot(edge3.cross(p.sub(p3))) < eps;
  
    if (side1 && side2 && side3) return new Intersection(ray, t, normal());

    return null;
  }
}

// Axis-aligned bounding box
class BBox extends Geometry {
  final Vec3 min, max;

  BBox(Vec3 min, Vec3 max) {
    this.min = min;
    this.max = max;
  }

  // Approximate normal based on the closest axis-aligned plane.
  Vec3 normal(Vec3 p) {
    final float eps = 1e-4;

    if (Math.abs(p.x - min.x) < eps) return new Vec3(-1, 0, 0);
    if (Math.abs(p.x - max.x) < eps) return new Vec3(1, 0, 0);
    if (Math.abs(p.y - min.y) < eps) return new Vec3(0, -1, 0);
    if (Math.abs(p.y - max.y) < eps) return new Vec3(0, 1, 0);
    if (Math.abs(p.z - min.z) < eps) return new Vec3(0, 0, 1);
    if (Math.abs(p.z - max.z) < eps) return new Vec3(0, 0, -1);

    return null;
  }

  /*
  * Based on PBRTv4's
  * [Bounds3::IntersectP](https://github.com/mmp/pbrt-v4/blob/39e01e61f8de07b99859df04b271a02a53d9aeb2/src/pbrt/util/vecmath.h#L1546C1-L1572C1)
  */
  Intersection intersect(Ray ray) {
    final float eps = 1e-5;

    float t0 = 0, t1 = Float.MAX_VALUE;
    for (int i = 0; i < 3; ++i) {
      final float d = ray.direction.at(i), o = ray.origin.at(i);
      final float min_val = min.at(i), max_val = max.at(i);
      final float invD = 1.f / d;
      float tNear = (min_val - o) * invD, tFar = (max_val - o) * invD;
      if (tNear > tFar) {
        float temp = tNear;
        tNear = tFar;
        tFar = temp;
      }

      tFar *= 1 + 2*eps;
      t0 = Math.max(t0, tNear);
      t1 = Math.min(t1, tFar);
      if (t0 > t1) return null;
    }

    final Vec3 norm = normal(ray.interp(t0));
    if (norm == null) throw new RuntimeException("Ray intersected BBox but could not determine the normal based on the intersection point.");

    return new Intersection(ray, t0, norm);
  }
}

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
  Hit(Intersection intersection, Surface surface) {
    this(intersection.t, intersection.point, intersection.normal, surface);
  }
}

// All objects rendered in a `Scene` inherit from `Object`.
abstract class Object {  
  abstract Hit raycast(Ray ray);
}

abstract class GeometryObject extends Object {
  final Geometry geometry;
  final Surface surface;

  GeometryObject(Geometry geometry, Surface surface) {
    this.geometry = geometry;
    this.surface = surface;
  }

  Hit raycast(Ray ray) {
    final Intersection intersection = geometry.intersect(ray);
    if (intersection != null) return new Hit(intersection, surface);
    return null;
  }
}

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
}

class TriangleObject extends GeometryObject {
  TriangleObject(Triangle triangle, Surface surface) {
    super(triangle, surface);
  }
}

class BBoxObject extends GeometryObject {
  BBoxObject(BBox bbox, Surface surface) {
    super(bbox, surface);
  }
}

// Bounding volume heirarchy acceleration data structure.
class BhvObject extends Object {
  class Node {
    final BBox bbox;
    final List<Object> objects;
    final List<Node> children;

    Node(BBox bbox, List<Object> objects) {
      this.bbox = bbox;
      this.objects = objects;
      this.children = new ArrayList();
    }

    Hit raycast(Ray ray) {
      Intersection intersection = bbox.intersect(ray);
      if (intersection == null) return null;
      for (Node c : children) {
        final Hit hit = c.raycast(ray);
        if (hit != null) return hit;
      }
      return null;
    }
  }

  Node root;
  Surface surface;

  BhvObject(BBox bbox, List<Object> objects, Surface surface) {
     // todo construct bvh
    root = new Node(bbox, objects);
    this.surface = surface;
  }

  Hit raycast(Ray ray) { return root.raycast(ray); }
}

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

  // Active accumulated triangle state:
  Surface surface = null;
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

  void clearVertices() { tri_a = tri_b = tri_c = null; }
  void commitVertices() {
    if (surface == null || tri_a == null || tri_b == null || tri_c == null) {
      println("Warning: Encountered an `End` command without a surface and three points. No triangle added.");
      return;
    }

    addObject(new TriangleObject(new Triangle(tri_a, tri_b, tri_c), surface));
  }

  void addBBox(BBox box) {
    final Mat4 transform = stack.top();
    addObject(new BBoxObject(new BBox(transform.transform(box.min), transform.transform(box.max)), surface));
  }

  void nameLatestObject(String name) {
    final Object object = popObject();
    if (object == null) {
      println("Warning: Attempted to name the latest added object " + name + ", but no " +
        (accelObjects != null ? "accellerated " : "") + "objects have been added to the scene.");
      return;
    }

    namedObjects.put(name, object);
  }

  // Create an instance of a named object and add that object to the list of scene objects.
  // Save the inverse of the current transformation matrix as part of the instance.
  void createInstance(String name) {
    if (!namedObjects.containsKey(name)) {
      println("Warning: Attempted to instance an object named " + name + ", but no object has been given this name.");
      return;
    }

    addObject(new InstancedObject(namedObjects.get(name), surface, name, stack.top()));
  }

  void beginAccel() { accelObjects = new ArrayList(); }
  void endAccel() { accelObjects = null; }

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
    final List<Object> activeObjects = accelObjects != null ? accelObjects : objects;
    activeObjects.add(object);
  }
  private Object popObject() {
    final List<Object> activeObjects = accelObjects != null ? accelObjects : objects;
    if (activeObjects.isEmpty()) return null;
    return activeObjects.remove(activeObjects.size() - 1);
  }
}

/**** Scene commands/parsing ****/

enum SceneCommandType {
  Background,
  Fov,
  Light,
  Surface,
  Begin,
  Vertex,
  End,

  Translate,
  Scale,
  Rotate,
  Push,
  Pop,

  Box,
  NamedObject,
  Instance,
  BeginAccel,
  EndAccel,

  Read,
  Render
}

// A `SceneCommand` is a complete, immutable, structured parsing of a textual scene command (a line in a `.cli` file).
// All commands with data provide a single `get()` method returning an instance of their respective data type.
abstract class SceneCommand {
  final SceneCommandType type;

  SceneCommand(SceneCommandType type) {
    this.type = type;
  }
}

class BackgroundCommand extends SceneCommand {
  final Color c;

  BackgroundCommand(float r, float g, float b) {
    super(SceneCommandType.Background);
    this.c = new Color(r, g, b);
  }

  Color get() { return c; }
}

class FovCommand extends SceneCommand {
  final float degrees;

  FovCommand(float degrees) {
    super( SceneCommandType.Fov);
    this.degrees = degrees;
  }

  float get() { return degrees; }
}

class LightCommand extends SceneCommand {
  final Light light;

  LightCommand(float x, float y, float z, float r, float g, float b) {
    super(SceneCommandType.Light);
    this.light = new Light(new Vec3(x, y, z), new Color(r, g, b));
  }

  Light get() { return light; }
}

class SurfaceCommand extends SceneCommand {
  final Surface surface;

  SurfaceCommand(float dr, float dg, float db) {
    super(SceneCommandType.Surface);
    this.surface = new Surface(new Color(dr, dg, db));
  }
  
  Surface get() { return surface; }
}

class BeginCommand extends SceneCommand {
  BeginCommand() { super(SceneCommandType.Begin); }
}

class VertexCommand extends SceneCommand {
  final Vec3 position;

  VertexCommand(float x, float y, float z) {
    super(SceneCommandType.Vertex);
    this.position = new Vec3(x, y, z);
  }

  Vec3 get() { return position; }
}

class EndCommand extends SceneCommand {
  EndCommand() { super(SceneCommandType.End); }
}

abstract class TransformCommand extends SceneCommand {
  final Mat4 transform;

  TransformCommand(SceneCommandType type, Mat4 transform) {
    super(type);
    this.transform = transform;
  }

  Mat4 get() { return transform; }
}

class TranslateCommand extends TransformCommand {
  TranslateCommand(float tx, float ty, float tz) {
    super(SceneCommandType.Translate, Mat4.translate(tx, ty, tz));
  }
}
class ScaleCommand extends TransformCommand {
  ScaleCommand(float sx, float sy, float sz) {
    super(SceneCommandType.Scale, Mat4.scale(sx, sy, sz));
  }
}
class RotateCommand extends TransformCommand {
  RotateCommand(float degrees, boolean x, boolean y, boolean z) {
    super(SceneCommandType.Rotate, Mat4.rotate(radians(degrees), x, y, z));
  }
}

class PushCommand extends SceneCommand {
  PushCommand() { super(SceneCommandType.Push); }
}
class PopCommand extends SceneCommand {
  PopCommand() { super(SceneCommandType.Pop); }
}

class BoxCommand extends SceneCommand {
  final BBox box;

  BoxCommand(float xmin, float ymin, float zmin, float xmax, float ymax, float zmax) {
    super(SceneCommandType.Box);
    this.box = new BBox(new Vec3(xmin, ymin, zmin), new Vec3(xmax, ymax, zmax));
  }
  
  BBox get() { return box; }
}

class NamedObjectCommand extends SceneCommand {
  final String name;

  NamedObjectCommand(String name) {
    super(SceneCommandType.NamedObject);
    this.name = name;
  }

  String get() { return name; }
}

class InstanceCommand extends SceneCommand {
  final String name;

  InstanceCommand(String name) {
    super(SceneCommandType.Instance);
    this.name = name;
  }

  String get() { return name; }
}

class BeginAccelCommand extends SceneCommand {
  BeginAccelCommand() { super(SceneCommandType.Begin); }
}
class EndAccelCommand extends SceneCommand {
  EndAccelCommand() { super(SceneCommandType.End); }
}

class ReadCommand extends SceneCommand {
  final String fileName;

  ReadCommand(String fileName) {
    super(SceneCommandType.Read);
    this.fileName = fileName;
  }

  String get() { return fileName; }
}

class RenderCommand extends SceneCommand {
  RenderCommand() { super(SceneCommandType.Render); }
}

class SceneCommandParser {
  SceneCommand parseTokens(String[] ts) {
    if (ts.length == 0) return null;

    final String name = ts.length == 0 ? "none" : ts[0];
    final SceneCommandType type = getType(name);
    if (type == null) {
      println("Warning: Unknown command type: " + name);
      return null;
    }

    switch (type) {
      case Background: return new BackgroundCommand(float(ts[1]), float(ts[2]), float(ts[3]));
      case Fov: return new FovCommand(float(ts[1]));
      case Light: return new LightCommand(float(ts[1]), float(ts[2]), float(ts[3]), float(ts[4]), float(ts[5]), float(ts[6]));
      case Surface: return new SurfaceCommand(float(ts[1]), float(ts[2]), float(ts[3]));
      case Begin: return new BeginCommand();
      case Vertex: return new VertexCommand(float(ts[1]), float(ts[2]), float(ts[3]));
      case End: return new EndCommand();

      case Translate: return new TranslateCommand(float(ts[1]), float(ts[2]), float(ts[3]));
      case Scale: return new ScaleCommand(float(ts[1]), float(ts[2]), float(ts[3]));
      case Rotate: return new RotateCommand(float(ts[1]), int(ts[2]) == 1, int(ts[3]) == 1, int(ts[4]) == 1);
      case Push: return new PushCommand();
      case Pop: return new PopCommand();

      case Box: return new BoxCommand(float(ts[1]), float(ts[2]), float(ts[3]), float(ts[4]), float(ts[5]), float(ts[6]));
      case NamedObject: return new NamedObjectCommand(ts[1]);
      case Instance: return new InstanceCommand(ts[1]);
      case BeginAccel: return new BeginAccelCommand();
      case EndAccel: return new EndAccelCommand();

      case Read: return new ReadCommand(ts[1]);
      case Render: return new RenderCommand();
      default: return null;
    }
  }

  SceneCommand parseTokens(String tokens) { return parseTokens(splitTokens(tokens, " ")); }

  // Assumes `filePath` is relative to `./data/`.
  List<SceneCommand> parseFile(String filePath) {
    final String[] lines = loadStrings(filePath);
    if (lines == null) {
      println("Error! Failed to read the file " + filePath);
      return new ArrayList();
    }

    return Arrays.stream(lines)
      .filter(line -> !line.trim().isEmpty() && !line.trim().startsWith("#")) // Filter out empty lines and comments.
      .map(this::parseTokens).filter(Objects::nonNull).collect(Collectors.toList());
  }

  // Returns null if there is no matching `SceneCommandType`.
  SceneCommandType getType(String name) {
    switch (name.toLowerCase()) {
      case "background": return SceneCommandType.Background;
      case "fov": return SceneCommandType.Fov;
      case "light": return SceneCommandType.Light;
      case "surface": return SceneCommandType.Surface;

      case "begin": return SceneCommandType.Begin;
      case "vertex": return SceneCommandType.Vertex;
      case "end": return SceneCommandType.End;

      case "translate": return SceneCommandType.Translate;
      case "scale": return SceneCommandType.Scale;
      case "rotate": return SceneCommandType.Rotate;
      case "push": return SceneCommandType.Push;
      case "pop": return SceneCommandType.Pop;

      case "box": return SceneCommandType.Box;
      case "named_object": return SceneCommandType.NamedObject;
      case "instance": return SceneCommandType.Instance;
      case "begin_accel": return SceneCommandType.BeginAccel;
      case "end_accel": return SceneCommandType.EndAccel;

      case "read": return SceneCommandType.Read;
      case "render": return SceneCommandType.Render;

      default: return null;
    }
  }
};

// Parse the text in a scene description file into a list of `SceneCommand`s.
// Then, iterate through the commands, updating and drawing the provided scene according to the parsed commands.
void interpret(String filePath, Scene scene) {
  final SceneCommandParser parser = new SceneCommandParser();
  final List<SceneCommand> commands = parser.parseFile(filePath);

  for (SceneCommand command : commands) {
    switch (command.type) {
      case Background:
        scene.backgroundColor = ((BackgroundCommand)command).get();
        break;
      case Fov:
        scene.fovDegrees = ((FovCommand)command).get();
        break;
      case Light:
        scene.lights.add(((LightCommand)command).get());
        break;
      case Surface:
        scene.surface = ((SurfaceCommand)command).get();
        break;
      case Begin:
        scene.clearVertices();
        break;
      case Vertex:
        scene.addVertex(((VertexCommand)command).get());
        break;
      case End:
        scene.commitVertices();
        break;
      case Translate:
      case Scale:
      case Rotate:
        scene.stack.apply(((TransformCommand)command).get());
        break;
      case Push:
        scene.stack.push();
        break;
      case Pop:
        scene.stack.pop();
        break;
      case Box:
        scene.addBBox(((BoxCommand)command).get());
        break;
      case NamedObject:
        scene.nameLatestObject(((NamedObjectCommand)command).get());
        break;
      case Instance:
        scene.createInstance(((InstanceCommand)command).get());
        break;
      case BeginAccel:
        scene.beginAccel();
        break;
      case EndAccel:
        scene.endAccel();
        break;
      case Read:
        interpret(((ReadCommand)command).get(), scene);
        break;
      case Render:
        drawScene(scene);
        break;
    }
  }
}

// Compute diffuse color at the hit point using the scene's light sources.
// Returns the scene's background coler if there is no hit.
Color shadeDiffuse(Hit hit, Scene scene) {
  if (hit == null) return scene.backgroundColor;

  final Vec3 N = hit.normal, P = hit.point;
  final Color diffuse = hit.surface.diffuse;
  return scene.lights.stream()
    .filter(light -> {
      // If the shadow ray intersects a scene object _before_ it hits the light,
      // the light does _not_ contribute to this point.
      final Vec3 pToL = light.position.sub(P), pToLDir = pToL.normalize();
      // Start the shadow ray epsilon away from the surface to prevent "immediately" hitting the surface _at_ `P`.
      final Ray shadowRay = new Ray(P.add(pToLDir.mult(1e-4)), pToLDir);
      final Hit shadowHit = scene.raycast(shadowRay);
      return shadowHit == null || shadowHit.t >= pToL.length();
    })
    .map(light -> {
      final Vec3 lightDir = light.position.sub(P).normalize();
      return diffuse.mult(light.c).mult(max(N.dot(lightDir), 0));
    })
    .reduce(new Color(0, 0, 0), Color::add);
}

void drawScene(Scene scene) {
  final Vec3 cameraPos = new Vec3(0, 0, 0);
  final float kw = tan(radians(scene.fovDegrees) / 2);
  final float kh = kw * (float(height) / float(width)); // Scale by aspect ratio.
  for(int y = 0; y < height; y++) {
    for(int x = 0; x < width; x++) {
      final Vec3 viewPlanePos = new Vec3(
        (x - width / 2.0) * (2 * kw / width),
        (y - height / 2.0) * (-2 * kh / height),
        -1
      );
      // This ray starts at the camera and points at the pixel on the view plane.
      final Ray cameraRay = new Ray(cameraPos, viewPlanePos.normalize());
      final Hit hit = scene.raycast(cameraRay);
      final Color c = shadeDiffuse(hit, scene);
      set(x, y, c.get()); // Set the color of the pixel
    }
  }
}

// prints mouse location clicks, for help in debugging
void mousePressed() {
  println ("You pressed the mouse at " + mouseX + " " + mouseY);
}

// you don't need to add anything in the "draw" function for this project
void draw() {}
