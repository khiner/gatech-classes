// Create polygonalized implicit surfaces.

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import java.lang.FunctionalInterface;

/** Specific implicit functions and helpers **/

// ImplicitInterface sphere = (x, y, z) -> sqrt(x*x + y*y + z*z);

class Ellipsoid extends SurfaceInstance {
  Ellipsoid(PVector center, PVector radii, PVector col) {
    super((x, y, z) -> {
      float dx = (center.x - x)/radii.x, dy = (center.y - y)/radii.y, dz = (center.z - z)/radii.z;
      return sqrt(dx*dx + dy*dy + dz*dz);
    }, col);
  }

  Ellipsoid(PVector center, float radius, PVector col) { this(center, new PVector(radius, radius, radius), col); }
}

class Vec2 {
  final float x, y;

  Vec2() { this.x = this.y = 0; }
  Vec2(float x, float y) {
    this.x = x;
    this.y = y;
  }

  Vec2 add(Vec2 o) { return new Vec2(x + o.x, y + o.y); }
  Vec2 sub(Vec2 o) { return new Vec2(x - o.x, y - o.y); }
  Vec2 mult(float scalar) { return new Vec2(x*scalar, y*scalar); }
  Vec2 div(float scalar) { return new Vec2(x/scalar, y/scalar); } 
  Vec2 flipX() { return new Vec2(-this.x, this.y); }
  Vec2 flipY() { return new Vec2(this.x, -this.y); }

  float length() { return sqrt(x*x + y*y); }
  Vec2 normalized() { return div(length()); }
}

Vec2 prev_mouse_pos = new Vec2();
PMatrix3D rot_mat;

Random rand = new Random();
PVector randVec(float scale) {
  return new PVector(
    rand.nextFloat()*scale - scale/2,
    rand.nextFloat()*scale - scale/2,
    rand.nextFloat()*scale - scale/2
  );
}

// Based on https://stackoverflow.com/a/7898685/780425
// (I hope it's fine to include this external code since it's not relevant to the subject matter!)
PVector hsvToRgb(float hue, float saturation, float v) {
    int h = (int)(hue*6);
    float f = hue*6 - h;
    float p = v*(1 - saturation);
    float q = v*(1 - f*saturation);
    float t = v*(1 - (1 - f)*saturation);
    switch (h) {
      case 0: return new PVector(v, t, p);
      case 1: return new PVector(q, v, p);
      case 2: return new PVector(p, v, t);
      case 3: return new PVector(p, q, v);
      case 4: return new PVector(t, p, v);
      case 5: return new PVector(v, p, q);
      default: throw new RuntimeException("Something went wrong when converting from HSV to RGB. Input was " + hue + ", " + saturation + ", " + v);
    }
}
PVector randCol() { return hsvToRgb(rand.nextFloat(), 0.4, 1.0); }

// camera parameters
final float CAMERA_DISTANCE_DEFAULT = 6;
float camera_distance = CAMERA_DISTANCE_DEFAULT;

class DrawFlags {
  boolean edges = false; // draw the polygon edges?
  boolean smooth_normals = false; // use smooth normals during shading?
}

DrawFlags draw_flags = new DrawFlags();
float iso_surface_threshold = 1;

void setup() {
  size(750, 750, OPENGL);
  
  // set up the rotation matrix
  rot_mat = (PMatrix3D)getMatrix();
  rot_mat.reset();
  
  // specify our implicit function is that of a sphere, then do isosurface extraction
  instances = List.of(new Ellipsoid(new PVector(0, 0, 0), 1, new PVector(1, 1, 1)));
  iso_surface_threshold = 1.0;
  isosurface();
}

void draw() {
  background(100, 100, 180); // clear the screen

  perspective(0.2*PI, 1.0, 0.01, 1000.0);
  camera(0, 0, camera_distance, 0, 0, 0, 0, 1, 0); // place the camera in the scene

  // create directional light sources
  directionalLight(100, 100, 100, -0.7, 0.7, -1);
  directionalLight(182, 182, 182, 0, 0, -1);
  
  pushMatrix();

  // decide if we are going to draw the polygon edges
  if (draw_flags.edges) stroke(0); // black edges
  else noStroke(); // no edges

  fill(250, 250, 250); // set the polygon color to white
  ambient(200, 200, 200);
  specular(0, 0, 0); // turn off specular highlights
  shininess(1.0);
  applyMatrix(rot_mat); // rotate the object using the global rotation matrix
  triangles.draw(); // draw the polygons from the implicit surface
  
  popMatrix();
}

// remember where the user clicked
void mousePressed() { prev_mouse_pos = new Vec2(mouseX, mouseY); }

// change the object rotation matrix while the mouse is being dragged
void mouseDragged() {
  if (!mousePressed) return;

  final Vec2 new_mouse_pos = new Vec2(mouseX, mouseY);
  final Vec2 mouse_delta = new_mouse_pos.sub(prev_mouse_pos).flipY();
  prev_mouse_pos = new_mouse_pos;
  if (mouse_delta.x == 0 && mouse_delta.y == 0) return;

  PMatrix3D rmat = new PMatrix3D();
  // Note: The rotation axes do not need to be normalized.
  rmat.rotate(0.005*mouse_delta.length(), mouse_delta.y, mouse_delta.x, 0);
  rot_mat.preApply(rmat);
}

// handle keystrokes
void keyPressed() {
  if (key == CODED) {
    if (keyCode == UP) camera_distance *= 0.9;
    else if (keyCode == DOWN) camera_distance /= 0.9;
    return;
  }

  switch(key) {
    case 'e': draw_flags.edges = !draw_flags.edges; break;
    case 'n': draw_flags.smooth_normals = !draw_flags.smooth_normals; break;
    case 'r': // reset camera view and rotation
      rot_mat.reset();
      camera_distance = CAMERA_DISTANCE_DEFAULT;
      break;
    case 'w': // write triangles to a file
      String filename = "implicit_mesh.cli";
      triangles.write(filename);
      println("wrote triangles to file: " + filename);
      break;
    case ',': // decrease the grid resolution
      if (gsize > 10) {
        gsize -= 10;
        isosurface();
      }
      break;
    case '.': // increase the grid resolution
      gsize += 10;
      isosurface();
      break;
    case '1':
      iso_surface_threshold = 1.0;
      instances = List.of(new Ellipsoid(new PVector(0, 0, 0), 1, new PVector(1, 1, 1)));
      isosurface();
      break;
    case '!':
      iso_surface_threshold = 1.0;
      camera_distance = CAMERA_DISTANCE_DEFAULT;
      rot_mat.reset();
      rot_mat.rotate(PI/6, -1, 0, 0);
      rot_mat.rotate(PI/18, 0, 0, -1);

      instances = List.of(new Ellipsoid(new PVector(0, 0, 0), new PVector(1, 1.0/3.0, 1), new PVector(1, 1, 1)));
      isosurface();
      break;
    case '2':
      iso_surface_threshold = 0.2;
      camera_distance = CAMERA_DISTANCE_DEFAULT;
      rot_mat.reset();

      // Draw pairs of blobby spheres that are partially blended together.
      // Draw one pair so that they are close together, and another where they have a small bridge connecting them.
      instances = Stream.of(
          new PVector(-0.6f, -1f, 0f), new PVector(0.6f, -1f, 0f),
          new PVector(-0.7f, 1f, 0f), new PVector(0.7f, 1f, 0f))
        .map(p -> new Ellipsoid(p, 1.0, new PVector(1, 1, 1)))
        .collect(Collectors.toList());
      isosurface();
      break;
    case '@':
      iso_surface_threshold = 0.25;

      // Draw 10 randomly placed/colored blobby spheres.
      instances = Stream.generate(() -> new Ellipsoid(randVec(3), 0.7, randCol())).limit(10).collect(Collectors.toList());
      isosurface();
      break;
  }
}

int timer_start_ms;
void resetTimer() { timer_start_ms = millis(); }
void printTimer() { println("elapsed ms: " + ((millis() - timer_start_ms) / 1000.0)); }
