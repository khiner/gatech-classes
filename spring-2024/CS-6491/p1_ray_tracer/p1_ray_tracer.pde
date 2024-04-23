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
      scene.lights.add(new PointLight(new Vec3(float(ts[1]), float(ts[2]), float(ts[3])), new Color(float(ts[4]), float(ts[5]), float(ts[6]))));
      break;
    // Create a disk light source, with center (x, y, z), a given radius, the direction that the light is facing (dx, dy, dz), and the light color (r, g, b).
    // Shadow rays are shot to random locations on this disk. When many rays per pixel are used, this creates soft shadows.
    case "disk_light":
      scene.lights.add(new DiskLight(float(ts[4]), new Vec3(float(ts[1]), float(ts[2]), float(ts[3])), new Vec3(float(ts[5]), float(ts[6]), float(ts[7])), new Color(float(ts[8]), float(ts[9]), float(ts[10]))));
      break;
    case "surface":
      scene.surface = new Surface(new Color(float(ts[1]), float(ts[2]), float(ts[3])));
      break;
    case "glossy":
      scene.surface = new Surface(new Color(float(ts[1]), float(ts[2]), float(ts[3])), new Color(float(ts[4]), float(ts[5]), float(ts[6])), float(ts[7]), float(ts[8]), float(ts[9]));
      break;
    case "rays_per_pixel":
      scene.raysPerPixel = int(ts[1]);
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
    // Change the last object that was defined into a moving object.
    // The values specify the amount of translation this object undergoes during one frame.
    // You can implement this in a manner similar to an instanced object, by translating the incoming ray's origin.
    // Remember that any shadow rays shot from this object should use the same random time!
    case "moving_object":
      scene.setLatestObjectVelocity(new Vec3(float(ts[1]), float(ts[2]), float(ts[3])));
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

// Compute color at the hit point using the scene's light sources and materials properties. //<>//
// Returns the scene's background color if there is no hit.
Color shade(Hit hit, Scene scene, float time) {
  if (hit == null) return scene.backgroundColor; //<>//

  final Vec3 N = hit.normal, P = hit.point;
  final Color diffuseColor = hit.surface.diffuse;
  Color c = new Color(0, 0, 0); // Accumulated return color
  for (Light light : scene.lights) {
    final Vec3 pToL = light.samplePosition().sub(P);
    final Vec3 pToLDir = pToL.normalize();
    // Start the shadow ray epsilo%5n away from the surface to prevent "immediately" hitting the surface at `P`.
    final Ray shadowRay = new Ray(P.add(pToLDir.mult(1e-4)), pToLDir, time);
    final Hit shadowHit = scene.raycast(shadowRay);
    if (shadowHit == null || shadowHit.t >= pToL.length()) {
      // Diffuse contribution
      final float diffuseIntensity = Math.max(N.dot(pToLDir), 0);
      c = c.add(diffuseColor.mult(light.c).mult(diffuseIntensity));
      // Specular contribution
      if (hit.surface.specular != null) {
        final Vec3 V = scene.cameraPosition.sub(P).normalize(); // Vector to the viewer
        final Vec3 H = pToLDir.add(V).normalize(); // Halfway vector
        final float specIntensity = (float)Math.pow(Math.max(N.dot(H), 0), hit.surface.specPower);
        c = c.add(hit.surface.specular.mult(light.c).mult(specIntensity));
      }
    }
  }

  // Reflective contribution
  if (hit.surface.reflectivity > 0) {
    Vec3 R = reflect(P.sub(scene.cameraPosition).normalize(), N); // Reflect direction
    if (hit.surface.glossRadius > 0) R = R.add(randomInSphere(hit.surface.glossRadius).normalize()); // Add fuzz factor

    final Ray reflectRay = new Ray(P.add(R.mult(1e-4)), R, time);
    final Hit reflectHit = scene.raycast(reflectRay);
    if (reflectHit != null) {
      final Color reflected = shade(reflectHit, scene, time); // Recursively reflected the ray
      c = c.add(reflected.mult(hit.surface.reflectivity));
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
    p = new Vec3((float)Math.random(), (float)Math.random(), (float)Math.random()).mult(2).sub(1);
  } while (p.dot(p) >= 1.0); // Ensure the point is within the unit sphere.

  return p.mult(radius);
}
 //<>//
void drawScene(Scene scene) {
  println("Drawing scene");
  reset_timer();
  final Vec3 cameraPos = scene.cameraPosition;
  final int raysPerPixel = scene.raysPerPixel;
  final boolean shouldSubsample = raysPerPixel > 1;
  final float kw = tan(radians(scene.fovDegrees) / 2);
  final float kh = kw * (float(height) / float(width)); // Scale by aspect ratio.
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      Color c = new Color(0, 0, 0);
      for (int r = 0; r < raysPerPixel; r++) {
        // Calculate a random offset within the pixel to create the sub-pixel ray.
        // This ray starts at the camera and points at the pixel on the view plane.
        final float offsetX = shouldSubsample ? (float)Math.random() - 0.5f : 0;
        final float offsetY = shouldSubsample ? (float)Math.random() - 0.5f : 0;
        final Vec3 viewPlanePos = new Vec3((x - width/2 + offsetX)*(2*kw/width), (y - height/2 + offsetY)*(-2*kh/height), -1);
        final float time = (float)Math.random(); // Time during the frame's duration
        final Ray cameraRay = new Ray(cameraPos, viewPlanePos.normalize(), time);
        final Hit hit = scene.raycast(cameraRay);
        c = c.add(shade(hit, scene, time));
      }
      c = c.div(raysPerPixel).mult(255); // Average the color across sub-pixel rays and scale range from [0,1] to [0,255].
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
