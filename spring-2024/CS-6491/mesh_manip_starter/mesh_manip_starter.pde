// Polygon mesh manipulation starter code.

// For object rotation by mouse
int mouseX_old = 0, mouseY_old = 0;
PMatrix3D rot_mat;
Mesh mesh;

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
  
  ambient(200, 200, 200);
  specular(0, 0, 0); // Turn off specular highlights
  shininess(1.0);
  
  applyMatrix(rot_mat); // Rotate the object using the global rotation matrix

  // THIS IS WHERE YOU SHOULD DRAW YOUR MESH

  if (mesh != null) mesh.draw();
    
  popMatrix();
}

void mousePressed() {
  mouseX_old = mouseX;
  mouseY_old = mouseY;
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

void readMesh(String filename) {
  mesh = loadMesh(filename);
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
  else if (mesh != null) {
    if (key == 'f') mesh.toggleSmoothShading();
    else if (key == 'w') mesh.toggleRandomColors();
    else if (key == 'e') mesh.toggleEdges();
    else if (key == 'v') mesh.toggleEdgeDebug();
    else if (key == 'n') mesh.moveDebugEdge(EdgeMove.Next);
    else if (key == 'p') mesh.moveDebugEdge(EdgeMove.Previous);
    else if (key == 'o') mesh.moveDebugEdge(EdgeMove.Opposite);
    else if (key == 's') mesh.moveDebugEdge(EdgeMove.Swing);
    else if (key == 'd') mesh = mesh.createDual();
    else if (key == 'g') mesh = mesh.subdivideMidpoint();
    else if (key == 'c') mesh = mesh.subdivideCatmullClark();
    else if (key == 'r') mesh.addRandomNoise(); // in-place
    else if (key == 'l') {
      for (int i = 0; i < 40; ++i) mesh.smoothLaplacian(0.6);
    }
    else if (key == 't') {
      for (int i = 0; i < 40; ++i) mesh.smoothTaubin(0.6307, -0.67315);
    }
  }
}
