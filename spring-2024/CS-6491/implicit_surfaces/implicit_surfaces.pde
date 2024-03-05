// Create polygonalized implicit surfaces.

// for object rotation by mouse
class Vec2 {
  float x = 0, y = 0;
  Vec2() {}
  Vec2(float x, float y) {
    this.x = x;
    this.y = y;
  }

  Vec2 add(Vec2 o) { return new Vec2(x + o.x, y + o.y); }
  Vec2 sub(Vec2 o) { return new Vec2(x - o.x, y - o.y); }
  Vec2 mult(float scalar) { return new Vec2(x*scalar, y*scalar); }
  Vec2 div(float scalar) { return new Vec2(x/scalar, y/scalar); } 

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
int timer;

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
  drawSurface(); // draw the polygons from the implicit surface
  
  popMatrix();
}

// remember where the user clicked
void mousePressed() { prev_mouse_pos = new Vec2(mouseX, mouseY); }

// change the object rotation matrix while the mouse is being dragged
void mouseDragged() {
  if (!mousePressed) return;

  final Vec2 new_mouse_pos = new Vec2(mouseX, mouseY);
  Vec2 mouse_delta = new_mouse_pos.sub(prev_mouse_pos);
  mouse_delta.y *= -1;
  if (mouse_delta.x == 0 && mouse_delta.y == 0) return;

  prev_mouse_pos = new_mouse_pos;

  float mouse_delta_len = mouse_delta.length();
  mouse_delta = mouse_delta.normalized();
  PMatrix3D rmat = (PMatrix3D)getMatrix();
  rmat.reset();
  rmat.rotate(0.005*mouse_delta_len, mouse_delta.y, mouse_delta.x, 0);
  rot_mat.preApply(rmat);
}

// handle keystrokes
void keyPressed()
{
  if (key == CODED) {
    if (keyCode == UP) camera_distance *= 0.9;
    else if (keyCode == DOWN) camera_distance /= 0.9;
    return;
  }
  
  if (key == 'e') draw_flags.edges = !draw_flags.edges;
  if (key == 'n') draw_flags.smooth_normals = !draw_flags.smooth_normals;
  if (key == 'r') {  // reset camera view and rotation
    rot_mat.reset();
    camera_distance = CAMERA_DISTANCE_DEFAULT;
  }
  if (key == 'w') { // write triangles to a file
    String filename = "implicit_mesh.cli";
    writeTriangles(filename);
    println ("wrote triangles to file: " + filename);
  }
  if (key == ',' && gsize > 10) { // decrease the grid resolution
    gsize -= 10;
    isosurface();
  }
  if (key == '.') { // increase the grid resolution
    gsize += 10;
    isosurface();
  }
  if (key == '1') {
    iso_surface_threshold = 1.0;
    implicit_func = a_sphere;
    isosurface();
  }
  if (key == '2') {}
}

void resetTimer() { timer = millis(); }
void printTimer() { println("elapsed ms: " + ((millis() - timer) / 1000.0)); }
