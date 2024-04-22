// This is the starter code for the CS 6491 Ray Tracing project.
//
// The most important part of the code is the interpreter, which will help
// you parse the scene description (.cli) files.

import java.util.Arrays;
import java.util.stream.Stream;

void setup() {
  size (300, 300);
  noStroke();
  background (0, 0, 0);
}

// Exclusing 's' prefix and extension
static String keyToFileName(char ch) {
  if (ch >= '1' && ch <= '9') return String.format("%02da", ch - '0');

  switch (ch) {
    case '!': return "01b";
    case '@': return "02b";
    case '#': return "03b";
    case '$': return "04b";
    case '%': return "05b";
    case '^': return "06b";
    case '&': return "07b";
    case '*': return "08b";
    case '(': return "09b";
    default: return null;
  }
}

void keyPressed() {
  final String file_name = keyToFileName(key);
  if (file_name == null) return;

  interpret(String.format("s%s.cli", file_name), new Scene());
}

int timer;
void reset_timer() { timer = millis(); }
void print_timer() { println("timer = " + (millis() - timer)/1000.0 + "s"); }

// Parse the text in a scene description file into its commands.
// Then, iterate through the commands, updating and drawing the provided scene according to the parsed commands.
// Assumes `filePath` is relative to `./data/`.
void interpret(String filePath, Scene scene) {
  final String[] lines = loadStrings(filePath);
  if (lines == null) throw new IllegalArgumentException("Failed to read the file " + filePath);

  Arrays.stream(lines)
    // Filter out empty lines and comments.
    .filter(line -> !line.trim().isEmpty() && !line.trim().startsWith("#"))
    .map(tokens -> splitTokens(tokens, " "))
    .forEach(ts -> {
    final String name = ts[0];
    switch (name) {
    case "background":
      scene.backgroundColor = new Color(float(ts[1]), float(ts[2]), float(ts[3]));
      break;
    case "fov":
      scene.fovDegrees = float(ts[1]);
      break;
    case "light":
      scene.lights.add(new Light(new Vec3(float(ts[1]), float(ts[2]), float(ts[3])), new Color(float(ts[4]), float(ts[5]), float(ts[6]))));
      break;
    case "surface":
      scene.surface = new Surface(new Color(float(ts[1]), float(ts[2]), float(ts[3])));
      break;
    case "glossy":
      scene.surface = new Surface(new Color(float(ts[1]), float(ts[2]), float(ts[3])), new Color(float(ts[4]), float(ts[5]), float(ts[6])), float(ts[7]), float(ts[8]), float(ts[9]));
      break;
    case "begin":
      scene.clearVertices();
      break;
    case "vertex":
      scene.addVertex(new Vec3(float(ts[1]), float(ts[2]), float(ts[3])));
      break;
    case "end":
      scene.commitVertices();
      break;

    case "translate":
      scene.stack.apply(Mat4.translate(float(ts[1]), float(ts[2]), float(ts[3])));
      break;
    case "scale":
      scene.stack.apply(Mat4.scale(float(ts[1]), float(ts[2]), float(ts[3])));
      break;
    case "rotate":
      scene.stack.apply(Mat4.rotate(radians(float(ts[1])), int(ts[2]) == 1, int(ts[3]) == 1, int(ts[4]) == 1));
      break;
    case "push":
      scene.stack.push();
      break;
    case "pop":
      scene.stack.pop();
      break;

    case "box":
      scene.addBBox(new BBox(new Vec3(float(ts[1]), float(ts[2]), float(ts[3])), new Vec3(float(ts[4]), float(ts[5]), float(ts[6]))));
      break;
    case "sphere":
      scene.addSphere(new Sphere(float(ts[1]), new Vec3(float(ts[2]), float(ts[3]), float(ts[4]))));
      break;
    case "named_object":
      scene.nameLatestObject(ts[1]);
      break;
    case "instance":
      scene.createInstance(ts[1]);
      break;
    case "begin_accel":
      scene.beginAccel();
      break;
    case "end_accel":
      scene.endAccel();
      break;

    case "read":
      interpret(ts[1], scene);
      break;
    case "render":
      drawScene(scene);
      break;
    }
  });
}

// Compute color at the hit point using the scene's light sources and materials properties.
// Returns the scene's background color if there is no hit.
Color shade(Hit hit, Scene scene) {
  if (hit == null) return scene.backgroundColor;

  final Vec3 N = hit.normal, P = hit.point;
  final Color diffuseColor = hit.surface.diffuse;
  Color c = new Color(0, 0, 0);

  for (Light light : scene.lights) {
    Vec3 pToL = light.position.sub(P);
    Vec3 pToLDir = pToL.normalize();
    // Start the shadow ray epsilon away from the surface to prevent "immediately" hitting the surface at `P`.
    Ray shadowRay = new Ray(P.add(pToLDir.mult(1e-4)), pToLDir);
    Hit shadowHit = scene.raycast(shadowRay);
    if (shadowHit == null || shadowHit.t >= pToL.length()) { //<>//
      // Diffuse contribution
      float diffuseIntensity = Math.max(N.dot(pToLDir), 0);
      Color diffuse = diffuseColor.mult(light.c).mult(diffuseIntensity); //<>//
      // Specular contribution
      if (hit.surface.specular != null) {
        Vec3 V = scene.cameraPosition.sub(P).normalize(); // Vector to the viewer
        Vec3 H = pToLDir.add(V).normalize(); // Halfway vector
        float specIntensity = (float) Math.pow(Math.max(N.dot(H), 0), hit.surface.specPower);
        Color specular = hit.surface.specular.mult(light.c).mult(specIntensity);
        diffuse = diffuse.add(specular);
      }
      c = c.add(diffuse);
    }
  }

  // Reflective contribution
  if (hit.surface.reflectivity > 0) {
    Vec3 R = reflect(P.sub(scene.cameraPosition).normalize(), N); // Reflect direction
    if (hit.surface.glossRadius > 0) {
      R = R.add(randomInSphere(hit.surface.glossRadius).normalize()); // Add fuzz factor
    }
    Ray reflectRay = new Ray(P.add(R.mult(1e-4)), R);
    Hit reflectHit = scene.raycast(reflectRay);
    if (reflectHit != null) {
      Color reflectedColor = shade(reflectHit, scene); // Recursive call for the reflected ray
      c = c.add(reflectedColor.mult(hit.surface.reflectivity));
    }
  }

  return c;
}

// Calculate reflection direction
Vec3 reflect(Vec3 I, Vec3 N) { return I.sub(N.mult(2 * I.dot(N))); }

// Generate a random point inside a sphere of a given radius, using rejection sampling.
Vec3 randomInSphere(float radius) {
    Vec3 p;
    do {
        p = new Vec3((float)Math.random(), (float)Math.random(), (float)Math.random()).mult(2.0f).sub(new Vec3(1, 1, 1));
    } while (p.dot(p) >= 1.0); // Ensure the point is within the unit sphere
    return p.mult(radius);
}
 //<>//
void drawScene(Scene scene) {
  println("Drawing scene");
  reset_timer();
  final Vec3 cameraPos = new Vec3(0, 0, 0);
  final float kw = tan(radians(scene.fovDegrees) / 2);
  final float kh = kw * (float(height) / float(width)); // Scale by aspect ratio.
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      final Vec3 viewPlanePos = new Vec3((x - width/2.0) * (2*kw / width), (y - height/2.0) * (-2*kh / height), -1);
      // This ray starts at the camera and points at the pixel on the view plane.
      final Ray cameraRay = new Ray(cameraPos, viewPlanePos.normalize());
      final Hit hit = scene.raycast(cameraRay);
      final Color c = shade(hit, scene).mult(255);
      set(x, y, color(c.r, c.g, c.b));
    }
  }
  print_timer();
}

// prints mouse location clicks, for help in debugging
void mousePressed() {
  println ("You pressed the mouse at " + mouseX + " " + mouseY);
}

// you don't need to add anything in the "draw" function for this project
void draw() {}
