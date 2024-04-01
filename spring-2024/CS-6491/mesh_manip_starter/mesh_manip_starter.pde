// Polygon mesh manipulation starter code.

// For object rotation by mouse
int mouseX_old = 0, mouseY_old = 0;
PMatrix3D rot_mat;

// Camera parameters
final float DEFAULT_CAMERA_DISTANCE = 6.0;
float camera_distance = DEFAULT_CAMERA_DISTANCE;

void setup() {
  size(800, 800, OPENGL);
  rot_mat = (PMatrix3D)getMatrix();
  rot_mat.reset();
}

void draw() {
  background(130, 130, 220); // Clear the screen to black

  perspective(PI*0.2, 1.0, 0.01, 1000.0);
  camera(0, 0, camera_distance, 0, 0, 0, 0, 1, 0); // Place the camera in the scene
  ambientLight(52, 52, 52); // Create an ambient light source

  // Create two directional light sources
  lightSpecular(0, 0, 0);
  directionalLight(150, 150, 150, -0.7, 0.7, -1);
  directionalLight(152, 152, 152, 0, 0, -1);
  
  pushMatrix();

  stroke(0); // Draw polygons with black edges
  fill(200, 200, 200); // Set the polygon color to white
  
  ambient(200, 200, 200);
  specular(0, 0, 0); // Turn off specular highlights
  shininess(1.0);
  
  applyMatrix(rot_mat); // Rotate the object using the global rotation matrix
  
  // THIS IS WHERE YOU SHOULD DRAW YOUR MESH

  beginShape();
  vertex(-1, 1, 0);
  vertex(1, 1, 0);
  vertex(0, -1, 0);
  endShape(CLOSE);
    
  popMatrix();
}

// Change the object rotation matrix while the mouse is being dragged
void mouseDragged() {
  if (!mousePressed) return;

  float dx = mouseX - mouseX_old;
  float dy = -(mouseY - mouseY_old);
  float len = sqrt(dx*dx + dy*dy);
  if (len == 0) len = 1;

  dx /= len;
  dy /= len;
  PMatrix3D rmat = (PMatrix3D)getMatrix();
  rmat.reset();
  rmat.rotate(len * 0.005, dy, dx, 0);
  rot_mat.preApply(rmat);

  mouseX_old = mouseX;
  mouseY_old = mouseY;
}

// handle keystrokes
void keyPressed() {
  if (key == CODED) {
    if (keyCode == UP) camera_distance *= 0.9; // zoom in
    else if (keyCode == DOWN) camera_distance /= 0.9; // zoom out
    return;
  }

  if (key == 'R') {
    rot_mat.reset();
    camera_distance = DEFAULT_CAMERA_DISTANCE;
  } else if (key == '1') readMesh("octa.ply");
  else if (key == '2') readMesh("cube.ply");
  else if (key == '3') readMesh("icos.ply");
  else if (key == '4') readMesh("dodeca.ply");
  else if (key == '5') readMesh("star.ply");
  else if (key == '6') readMesh("torus.ply");
  else if (key == '7') readMesh("s.ply");
}
