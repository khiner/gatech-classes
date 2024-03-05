// Create polygonalized implicit surfaces.

import java.lang.FunctionalInterface;

// This is a functional interface for defining implicit functions.
@FunctionalInterface
interface ImplicitInterface {
  float at(float x, float y, float z);
}

// Implicit function for a sphere at the origin.
ImplicitInterface a_sphere = (x, y, z) -> sqrt(x*x + y*y + z*z);

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
  implicit_func = a_sphere;
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

  final Vec2 mouse_delta_norm = mouse_delta.normalized();
  PMatrix3D rmat = (PMatrix3D)getMatrix();
  rmat.reset();
  rmat.rotate(0.005*mouse_delta.length(), mouse_delta_norm.y, mouse_delta_norm.x, 0);
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
      implicit_func = a_sphere;
      isosurface();
      break;
    case '2': break;
  }
}

int timer_start_ms;
void resetTimer() { timer_start_ms = millis(); }
void printTimer() { println("elapsed ms: " + ((millis() - timer_start_ms) / 1000.0)); }
