// This is the starter code for the CS 6491 Ray Tracing project.
//
// The most important part of the code is the interpreter, which will help
// you parse the scene description (.cli) files.

import java.util.Arrays;
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

color unitColor(float r, float g, float b) { return color(r * 255, g * 255, b * 255); }

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

class Light {
  PVector position;
  color c;

  Light(PVector position, color c) {
    this.position = position;
    this.c = c;
  }
  Light() {
    this(new PVector(0, 0, 0), unitColor(1, 1, 1));
  }
}

class Surface {
  color diffuse;
  
  Surface(color diffuse) {
    this.diffuse = diffuse;
  }
  Surface() {
    this(color(0, 0, 0));
  }
}

abstract class SceneCommand {
  String name;
  SceneCommandType type;

  SceneCommand(String name, SceneCommandType type) {
    this.name = name;
    this.type = type;
  }
}

class BackgroundCommand extends SceneCommand {
  color c;

  BackgroundCommand(String name, SceneCommandType type, float r, float g, float b) {
    super(name, type);
    this.c = unitColor(r, g, b);
  }
  
  color get() { return c; }
}

class FovCommand extends SceneCommand {
  float degrees;

  FovCommand(String name, SceneCommandType type, float degrees) {
    super(name, type);
    this.degrees = degrees;
  }
  
  float get() { return degrees; }
}

class LightCommand extends SceneCommand {
  Light light;

  LightCommand(String name, SceneCommandType type, float x, float y, float z, float r, float g, float b) {
    super(name, type);
    this.light = new Light(new PVector(x, y, z), unitColor(r, g, b));
  }
  
  Light get() { return light; }
}

class SurfaceCommand extends SceneCommand {
  Surface surface;

  SurfaceCommand(String name, SceneCommandType type, float dr, float dg, float db) {
    super(name, type);
    this.surface = new Surface(unitColor(dr, dg, db));
  }
  
  Surface get() { return surface; }
}

class BeginCommand extends SceneCommand {
  BeginCommand(String name, SceneCommandType type) { super(name, type); }
}

class VertexCommand extends SceneCommand {
  PVector position;

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


  List<SceneCommand> parseFile(String filePath) {
    println("Parsing '" + filePath + "'");

    String[] lines = loadStrings(filePath);
    if (lines == null) {
      println("Error! Failed to read the file.");
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

class Triangle {
    PVector p1, p2, p3;

    Triangle(PVector p1, PVector p2, PVector p3) {
        this.p1 = p1;
        this.p2 = p2;
        this.p3 = p3;
    }
}

class Scene {
  float fovDegrees;
  color backgroundColor;
  Light light;
  Surface surface;
  List<Triangle> triangles;

  Scene() {
    fovDegrees = 0;
    backgroundColor = color(0, 0, 0);
    light = new Light();
    surface = new Surface();
    triangles = new ArrayList();
  }

  // Returns the closest intersecting triangle to the ray's origin,
  // or `null` if the ray does not intersect any triangles.
  Triangle raycastTriangle(Ray ray) {
    Triangle result = null;
    float closestT = Float.MAX_VALUE;
    for (Triangle triangle : triangles) {
      float t = rayTriangleIntersection(ray, triangle);
      if (t > 0 && t < closestT) {
        closestT = t;
        result = triangle;
      }
    }

    return result;
  }
}

// This routine parses the text in a scene description file into a list of `SceneCommand`s,
// creates a new scene and iterates through the commands, updating and drawing the scene.
void interpreter(String filePath) {
  SceneCommandParser parser = new SceneCommandParser();
  final List<SceneCommand> commands = parser.parseFile(filePath);
  Scene scene = new Scene();
  
  PVector tri_1 = null, tri_2 = null, tri_3 = null;

  for (SceneCommand command : commands) {
    switch (command.type) {
      case Background:
        scene.backgroundColor = ((BackgroundCommand)command).get();
        break;
      case Fov:
        scene.fovDegrees = ((FovCommand)command).get();
        break;
      case Light:
        scene.light = ((LightCommand)command).get();
        break;
      case Surface:
        scene.surface = ((SurfaceCommand)command).get();
        break;
      case Begin:
        tri_1 = tri_2 = tri_3 = null;
        break;
      case Vertex:
        VertexCommand vertexCommand = (VertexCommand)command;
          if (tri_1 == null) tri_1 = vertexCommand.get();
          else if (tri_2 == null) tri_2 = vertexCommand.get();
          else if (tri_3 == null) tri_3 = vertexCommand.get();
          else throw new IllegalArgumentException("More than three vertices within a single begin/end block.");
        break;
      case End:
        if (tri_1 != null && tri_2 != null && tri_3 != null) {
          scene.triangles.add(new Triangle(tri_1, tri_2, tri_3));
        }
        break;
      case Render:
        drawScene(scene);
        break;
    }
  }
}

class Ray {
  PVector origin, direction;
  
  Ray(PVector origin, PVector direction) {
    this.origin = origin;
    this.direction = direction;
  }
}

// Moller-Trumbore ray-triangle intersection.
float rayTriangleIntersection(Ray ray, Triangle tri) {
  final PVector edge1 = PVector.sub(tri.p2, tri.p1), edge2 = PVector.sub(tri.p3, tri.p1);
  final PVector h = ray.direction.cross(edge2);
  final float a = edge1.dot(h);
  if (a > -0.00001 && a < 0.00001) return -1; // Ray is parallel to the triangle.

  final float f = 1.0 / a;
  final PVector s = PVector.sub(ray.origin, tri.p1);
  final float u = f * (s.dot(h));
  if (u < 0.0 || u > 1.0) return -1;

  final PVector q = s.cross(edge1);
  final float v = f * ray.direction.dot(q);
  if (v < 0.0 || u + v > 1.0) return -1;

  // Compute t to find out where the intersection point is on the line.
  final float t = f * edge2.dot(q);
  if (t > 0.00001) return t; // Ray intersects.

  return -1; // Line intersects, but ray does not.
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
      final Triangle intersectingTriangle = scene.raycastTriangle(cameraRay);
      final color c = intersectingTriangle == null ? scene.backgroundColor : scene.surface.diffuse;
      // todo Compute diffuse color at the intersection point using the provided point light sources.
      /* Notes:      
        https://gatech.instructure.com/courses/360970/external_tools/18649
        One way to make sure the normals always facing outward is to make the normal and ray directions consistent. You can do something like:
        if dot(N,ray_dir) < 0) // Keep the same normal
        else // Flip normal direction
      */

      set(x, y, c); // Set the color of the pixel
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
