// Create polygonalized implicit surfaces.

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import java.lang.FunctionalInterface;

// Positioning, scaling, and rotation happens in `Primitive`.

ImplicitInterface sphere = p -> p.mag();

class LineSegment implements ImplicitInterface {
  final PVector a, b;

  LineSegment(PVector a, PVector b) {
    this.a = a;
    this.b = b;
  }

  float at(PVector p) {
    final PVector ap = PVector.sub(p, a), ab = PVector.sub(b, a);
    final float t = ap.dot(ab) / ab.magSq();
    final PVector closest = t < 0 ? a : t > 1 ? b : PVector.add(a, PVector.mult(ab, t));
    return PVector.sub(p, closest).mag();
  }
}

class Torus implements ImplicitInterface {
  final float R, r;

  Torus(float R, float r) {
    this.R = R;
    this.r = r;
  }

  float at(PVector p) {
    return (float)Math.pow(p.magSq() + R*R - r*r, 2) - 4*R*R*(vec(p.x, p.y).magSq());
  }
}

// Use x coord as twist angle.
class TwistX implements WarpInterface {
  final float k; // Frequency multiple (frequency is `kx`)

  TwistX(float k) { this.k = k; }

  PVector apply(PVector p) {
    return vec(p.x, p.y*cos(k*p.x) - p.z*sin(k*p.x), p.y*sin(k*p.x) + p.z*cos(k*p.x));
  }
}

// Scale based on x position.
class TaperX implements WarpInterface {
  final float k_1, k_2, x_min, x_max; // Tapering factor

  TaperX(float x_min, float x_max, float k_1, float k_2) {
    this.k_1 = k_1;
    this.k_2 = k_2;
    this.x_min = x_min;
    this.x_max = x_max;
  }

  PVector apply(PVector p) { return vec(p.x, p.y/k(p.x), p.z/k(p.x)); }

  private float t(float x) {
    if (x < x_min) return 0;
    if (x_min <= x  || x <= x_max) return (x - x_min)/(x_max - x_min);
    return 1;
  }

  private float k(float x) { return (1 - t(x))*k_1 + t(x)*k_2; }
}

abstract class BooleanCombination implements ImplicitInterface {
  ImplicitInterface a, b;

  BooleanCombination(ImplicitInterface a, ImplicitInterface b) {
    this.a = a;
    this.b = b;
  }
}

class ImplicitUnion extends ImplicitBoolean {
  ImplicitUnion(ImplicitInterface a, ImplicitInterface b) { super(a, b); }
  float at(PVector p) { return min(a.at(p), b.at(p)); }
}

class ImplicitIntersection extends ImplicitBoolean {
  ImplicitIntersection(ImplicitInterface a, ImplicitInterface b) { super(a, b); }
  float at(PVector p) { return max(a.at(p), b.at(p)); }
}

class ImplicitSubtraction extends ImplicitBoolean {
  ImplicitSubtraction(ImplicitInterface a, ImplicitInterface b) { super(a, b); }
  float at(PVector p) { return max(a.at(p), vmax-b.at(p)); }
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
}

Vec2 prev_mouse_pos = new Vec2();
PMatrix3D rot_mat;

Random rand = new Random();
PVector randVec(float scale) {
  return vec(
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
      case 0: return vec(v, t, p);
      case 1: return vec(q, v, p);
      case 2: return vec(p, v, t);
      case 3: return vec(p, q, v);
      case 4: return vec(t, p, v);
      case 5: return vec(v, p, q);
      default: throw new RuntimeException("Something went wrong when converting from HSV to RGB. Input was " + hue + ", " + saturation + ", " + v);
    }
}
PVector randCol() { return hsvToRgb(rand.nextFloat(), 0.4, 1.0); }

// Note: The rotation axes do _not_ need to be normalized.
PMatrix3D rotationMatrix(float angle, float x, float y, float z) {
  PMatrix3D rmat = new PMatrix3D();
  rmat.rotate(angle, x, y, z);
  return rmat;
}

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
  primitives = List.of(new Primitive(sphere));
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

  rot_mat.preApply(rotationMatrix(0.005*mouse_delta.length(), mouse_delta.y, mouse_delta.x, 0));
}

void resetScene(float threshold) {
  iso_surface_threshold = threshold;
  rot_mat.reset();
  camera_distance = CAMERA_DISTANCE_DEFAULT;
}
void resetScene() { resetScene(iso_surface_threshold); }

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
      resetScene();
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
      setPrimitive(new Primitive(sphere));
      break;
    case '!':
      resetScene(1.0);
      rot_mat.rotate(PI/6, -1, 0, 0);
      rot_mat.rotate(PI/18, 0, 0, -1);
      setPrimitive(new Primitive(sphere, DEFAULT_COL, vec(0, 0, 0), vec(1, 1.0/3.0, 1)));
      break;
    case '2':
      resetScene(0.2);
      // Draw pairs of blobby spheres that are partially blended together.
      // Draw one pair so that they are close together, and another where they have a small bridge connecting them.
      setPrimitives(Stream.of(
          vec(-0.6, -1, 0), vec(0.6, -1, 0),
          vec(-0.7, 1, 0), vec(0.7, 1, 0))
        .map(p -> new Primitive(sphere, DEFAULT_COL, p))
        .collect(Collectors.toList()));
      break;
    case '@':
      resetScene(0.4);
      // Draw 10 randomly placed/colored blobby spheres.
      setPrimitives(Stream.generate(() -> new Primitive(sphere, randCol(), randVec(3))).limit(10).collect(Collectors.toList()));
      break;
    case '3':
      resetScene(0.7);
      setPrimitive(new Primitive(new LineSegment(vec(-1, 0, 0), vec(1, 0, 0))));
      break;
    case '#':
      resetScene(0.4);
      camera_distance = CAMERA_DISTANCE_DEFAULT*1.5;
      setPrimitives(List.of(
        new Primitive(new LineSegment(vec(-1.2, -1.2, 0), vec(1.2, -1.2, 0))),
        new Primitive(new LineSegment(vec(1.2, -1.2, 0), vec(1.2, 1.2, 0))),
        new Primitive(new LineSegment(vec(1.2, 1.2, 0), vec(-1.2, 1.2, 0))),
        new Primitive(new LineSegment(vec(-1.2, 1.2, 0), vec(-1.2, -1.2, 0)))
      ));
      break;
    case '4':
      resetScene(0.1);
      setPrimitive(new Primitive(new Torus(1, 0.5)));
      break;
    case '$':
      resetScene(0.4);
      setPrimitives(List.of(
        new Primitive(new Torus(1.2, 0.0), DEFAULT_COL, vec(-1, 0, 0), 0.45, new Rotation(PI/4, 1, 0, 0)),
        new Primitive(new Torus(1.2, 0.0), DEFAULT_COL, vec(0, 0, 0), 0.45),
        new Primitive(new Torus(1.2, 0.0), DEFAULT_COL, vec(1, 0, 0), 0.45, new Rotation(-PI/4, 1, 0, 0))
      ));
      break;
    case '5':
      resetScene(0.3);
      setPrimitive(new Primitive(new LineSegment(vec(-2.3, 0.4, 0), vec(2.3, 0.4, 0)), 0.7));
      break;
    case '%':
      resetScene(0.3);
      setPrimitive(new Primitive(new LineSegment(vec(-2.3, 0.4, 0), vec(2.3, 0.4, 0)), 0.7).addWarp(new TwistX(5)));
      break;
    case '6':
      resetScene(0.3);
      setPrimitive(
        new Primitive(new LineSegment(vec(-1.2, 0.0, 0), vec(1.2, 0.0, 0)))
          .addWarp(new TaperX(-1.3, 1.3, 0.3, 1.1))
      );
      break;
    case '^':
      resetScene(0.3);
      setPrimitive(
        new Primitive(new LineSegment(vec(-1.2, 0.2, 0), vec(1.2, 0.2, 0)))
          .addWarp(new TwistX(10))
          .addWarp(new TaperX(-1.3, 1.3, 0.3, 1.1))
      );
      break;
    case '7':
      resetScene(1.0);
      rot_mat.rotate(PI/16, -1, 0, 0);
      rot_mat.rotate(PI/16, 0, 0, -1);
      setPrimitive(new Primitive(
        new ImplicitIntersection(
          new Primitive(sphere, DEFAULT_COL, vec(0, 0.5, 0)), new Primitive(sphere, DEFAULT_COL, vec(0, -0.5, 0))
        )
      ));
      break;
    case '&':
      resetScene(1.0);
      setPrimitive(new Primitive(
        new ImplicitSubtraction(
          sphere, new Primitive(new LineSegment(vec(0, 0, -vmax), vec(0, 0.0, vmax)), 0.5)
        )
      ));
      break;
    case '8':
      resetScene(1.0);
      setPrimitive(new Primitive(
        new ImplicitSubtraction(
          sphere, new Primitive(new LineSegment(vec(0, 0, -vmax), vec(0, 0.0, vmax)), 0.5)
        )
      ));
      break;
  }
}

int timer_start_ms;
void resetTimer() { timer_start_ms = millis(); }
void printTimer() { println("elapsed ms: " + ((millis() - timer_start_ms) / 1000.0)); }
