// This is the starter code for the CS 6491 Ray Tracing project.
//
// The most important part of the code is the interpreter, which will help
// you parse the scene description (.cli) files.

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

boolean debug_flag = false;

void setup() {
  size (300, 300);  
  noStroke();
  background (0, 0, 0);
}

void keyPressed() {
  switch(key) {
    case '1': interpreter("s1.cli"); break;
    case '2': interpreter("s2.cli"); break;
    case '3': interpreter("s3.cli"); break;
    case '4': interpreter("s4.cli"); break;
  }
}

/**** Scene description ****/

class Color {
  float r, g, b; // In range [0, 1]

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
  PVector position;
  Color c;

  Light(PVector position, Color c) {
    this.position = position;
    this.c = c;
  }
}

class Surface {
  Color diffuse;
  
  Surface(Color diffuse) {
    this.diffuse = diffuse;
  }
}

class Triangle {
  Surface surface;
  PVector p1, p2, p3;

  Triangle(Surface surface, PVector p1, PVector p2, PVector p3) {
    this.p1 = p1;
    this.p2 = p2;
    this.p3 = p3;
    this.surface = surface;
  }

  PVector normal() {
    final PVector edge1 = PVector.sub(p2, p1), edge2 = PVector.sub(p3, p1);
    return edge1.cross(edge2).normalize();
  }
}

class Ray {
  PVector origin, direction;
  
  Ray(PVector origin, PVector direction) {
    this.origin = origin;
    this.direction = direction;
  }
 
  // Interpolate along the ray to find the point at distance `t`.
  PVector interp(float t) {
    return PVector.add(origin, PVector.mult(direction, t));
  }
}

class Hit {
  Ray ray;
  float t; // Distance along the ray
  PVector point; // Intersection point
  Triangle triangle;
  
  Hit(Ray ray, float t, Triangle triangle) {
    this.ray = ray;
    this.t = t;
    this.triangle = triangle;
    this.point = ray.interp(t);
  }
}

class Scene {
  float fovDegrees;
  Color backgroundColor;
  List<Light> lights;
  List<Triangle> triangles;

  Scene() {
    fovDegrees = 0;
    backgroundColor = new Color(0, 0, 0);
    lights = new ArrayList();
    triangles = new ArrayList();
  }

  // Returns a hit for the intersecting triangle closest to the ray's origin,
  // or `null` if the ray does not intersect any triangle.
  Hit raycast(Ray ray) {
    return triangles.stream()
      .map(triangle -> rayTriangleIntersection(ray, triangle))
      .filter(Objects::nonNull)
      .min(Comparator.comparingDouble(hit -> hit.t))
      .orElse(null);
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
  Render
}

// A `SceneCommand` is a complete, immutable, structured parsing of a textual scene command (a line in a `.cli` file).
// All commands with data provide a single `get()` method returning an instance of their respective data type.
abstract class SceneCommand {
  final String name;
  final SceneCommandType type;

  SceneCommand(String name, SceneCommandType type) {
    this.name = name;
    this.type = type;
  }
}

class BackgroundCommand extends SceneCommand {
  final Color c;

  BackgroundCommand(String name, SceneCommandType type, float r, float g, float b) {
    super(name, type);
    this.c = new Color(r, g, b);
  }
  
  Color get() { return c; }
}

class FovCommand extends SceneCommand {
  final float degrees;

  FovCommand(String name, SceneCommandType type, float degrees) {
    super(name, type);
    this.degrees = degrees;
  }
  
  float get() { return degrees; }
}

class LightCommand extends SceneCommand {
  final Light light;

  LightCommand(String name, SceneCommandType type, float x, float y, float z, float r, float g, float b) {
    super(name, type);
    this.light = new Light(new PVector(x, y, z), new Color(r, g, b));
  }
  
  Light get() { return light; }
}

class SurfaceCommand extends SceneCommand {
  final Surface surface;

  SurfaceCommand(String name, SceneCommandType type, float dr, float dg, float db) {
    super(name, type);
    this.surface = new Surface(new Color(dr, dg, db));
  }
  
  Surface get() { return surface; }
}

class BeginCommand extends SceneCommand {
  BeginCommand(String name, SceneCommandType type) { super(name, type); }
}

class VertexCommand extends SceneCommand {
  final PVector position;

  VertexCommand(String name, SceneCommandType type, float x, float y, float z) {
    super(name, type);
    this.position = new PVector(x, y, z);
  }
  
  PVector get() { return position; }
}

class EndCommand extends SceneCommand {
  EndCommand(String name, SceneCommandType type) { super(name, type); }
}

class RenderCommand extends SceneCommand {
  RenderCommand(String name, SceneCommandType type) { super(name, type); }
}

class SceneCommandParser {
  SceneCommand parseTokens(String[] tokens) {
    if (tokens.length == 0) return null;

    final String name = tokens.length == 0 ? "none" : tokens[0];
    final SceneCommandType type = getType(name);
    switch (type) {
      case Background:
        return new BackgroundCommand(name, type, float(tokens[1]), float(tokens[2]), float(tokens[3]));
      case Fov:
        return new FovCommand(name, type, float(tokens[1]));
      case Light:
        return new LightCommand(name, type, float(tokens[1]), float(tokens[2]), float(tokens[3]), float(tokens[4]), float(tokens[5]), float(tokens[6]));
      case Surface:
        return new SurfaceCommand(name, type, float(tokens[1]), float(tokens[2]), float(tokens[3]));
      case Begin:
        return new BeginCommand(name, type);
      case Vertex:
        return new VertexCommand(name, type, float(tokens[1]), float(tokens[2]), float(tokens[3]));
      case End:
        return new EndCommand(name, type);
      case Render:
        return new RenderCommand(name, type);
      default:
        throw new IllegalArgumentException("Unknown command type: " + name);
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

  SceneCommandType getType(String name) {
    switch (name.toLowerCase()) {
      case "background": return SceneCommandType.Background;
      case "fov": return SceneCommandType.Fov;
      case "light": return SceneCommandType.Light;
      case "surface": return SceneCommandType.Surface;
      case "begin": return SceneCommandType.Begin;
      case "vertex": return SceneCommandType.Vertex;
      case "end": return SceneCommandType.End;
      case "render": return SceneCommandType.Render;
      default: throw new IllegalArgumentException("Unknown command type: " + name);
    }
  }
};

// This routine parses the text in a scene description file into a list of `SceneCommand`s,
// creates a new scene and iterates through the commands, updating and drawing the scene.
void interpreter(String filePath) {
  final SceneCommandParser parser = new SceneCommandParser();
  final List<SceneCommand> commands = parser.parseFile(filePath);

  // Mutable scene, populated and rendered according to the parsed commands.
  Scene scene = new Scene();
  // Active accumulated triangle state:
  Surface surface = null;
  PVector tri_a = null, tri_b = null, tri_c = null;

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
        surface = ((SurfaceCommand)command).get();
        break;
      case Begin:
        tri_a = tri_b = tri_c = null;
        break;
      case Vertex:
        final PVector vertex = ((VertexCommand)command).get();
        if (tri_a == null) tri_a = vertex;
        else if (tri_b == null) tri_b = vertex;
        else if (tri_c == null) tri_c = vertex;
        else throw new IllegalArgumentException("More than three vertices within a single begin/end block.");
        break;
      case End:
        if (surface != null && tri_a != null && tri_b != null && tri_c != null) {
          scene.triangles.add(new Triangle(surface, tri_a, tri_b, tri_c));
        } else {
          println("Warning: Encountered an `End` command without a surface and three points. No triangle added.");
        }
        break;
      case Render:
        drawScene(scene);
        break;
    }
  }
}

// Returns a hit representing the intersection, or `null` if there is no intersection.
Hit rayTriangleIntersection(Ray ray, Triangle tri) {
  final float eps = 0.00001;

  final PVector N = tri.normal();
  // Calculate `t` (the distance from the ray origin to the intersection point).
  final float denom = N.dot(ray.direction);
  if (abs(denom) < eps) return null; // Ray is parallel to the triangle.

  final float t = -(N.dot(ray.origin) - N.dot(tri.p1)) / denom;
  if (t < 0) return null; // Intersects behind the ray's origin.

  // Check if the intersection point is inside the triangle.
  final PVector p = ray.interp(t);
  final PVector edge1 = PVector.sub(tri.p2, tri.p1), edge2 = PVector.sub(tri.p3, tri.p2), edge3 = PVector.sub(tri.p1, tri.p3);
  final boolean
    side1 = -N.dot(edge1.cross(PVector.sub(p, tri.p1))) < eps,
    side2 = -N.dot(edge2.cross(PVector.sub(p, tri.p2))) < eps,
    side3 = -N.dot(edge3.cross(PVector.sub(p, tri.p3))) < eps;

  return side1 && side2 && side3 ? new Hit(ray, t, tri) : null;
}

// Compute diffuse color at the hit point using the scene's light sources.
// Returns the scene's background coler if there is no hit.
Color shadeDiffuse(Hit hit, Scene scene) {
  if (hit == null) return scene.backgroundColor;

  PVector N = hit.triangle.normal();
  if (N.dot(hit.ray.direction) > 0) N.mult(-1); // Ensure the normal is facing the camera.

  final Color diffuse = hit.triangle.surface.diffuse;
  return scene.lights.stream().map(light -> {
    final PVector lightDir = PVector.sub(light.position, hit.point).normalize();
    return diffuse.mult(light.c).mult(max(N.dot(lightDir), 0));
  }).reduce(new Color(0, 0, 0), Color::add);
}

void drawScene(Scene scene) {
  final PVector cameraPos = new PVector(0, 0, 0);
  final float kw = tan(radians(scene.fovDegrees) / 2);
  final float kh = kw * (float(height) / float(width)); // Scale by aspect ratio.
  for(int y = 0; y < height; y++) {
    for(int x = 0; x < width; x++) {
      final PVector viewPlanePos = new PVector(
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
void draw() {
}
