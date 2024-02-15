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

static Integer charToCliFileNumber(char ch) {
  if (ch >= '1' && ch <= '9') return ch - '0';
  if (ch == '0') return 10;
  if (ch >= 'a' && ch <= 'd') return ch - 'a' + 11;
  return null;
}

void keyPressed() {
  final Integer cli_file_number = charToCliFileNumber(key);
  if (cli_file_number == null) return;

  interpret(String.format("s%02d.cli", cli_file_number), new Scene());
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
      final Color c = shadeDiffuse(hit, scene).mult(255);
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
