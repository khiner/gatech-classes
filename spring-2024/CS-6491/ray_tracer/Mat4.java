import processing.core.PMatrix3D;

/**
4x4 matrix
Instances are immutable.
Rotation angles are in radians.
*/
public class Mat4 {
  private final float[] data; // Row-major, using 1D array instead of 2D for better cache locality.

  public Mat4() {
    data = new float[16];
    // Initialize identity matrix.
    data[0] = data[5] = data[10] = data[15] = 1;
  }

  public Mat4(PMatrix3D other) {
    data = new float[16];
    other.get(data);
  }

  // Copy constructor
  public Mat4(Mat4 other) { this.data = other.data.clone(); }

  public Mat4 invert() {
    PMatrix3D pmat4 = new PMatrix3D();
    pmat4.set(data);
    if (pmat4.invert()) return new Mat4(pmat4);
    return null; // Inversion not successful. (Some matrices map more than one point to the same image point, and so are irreversible.)
  }

  public Mat4 transpose() {
    Mat4 result = new Mat4();
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        result.data[i * 4 + j] = this.data[j * 4 + i];
      }
    }
    return result;
  }

  public Mat4 mult(Mat4 other) {
    Mat4 result = new Mat4();
    for (int row = 0; row < 4; row++) {
      for (int col = 0; col < 4; col++) {
        result.data[row * 4 + col] = 0;
        for (int k = 0; k < 4; k++) {
          result.data[row * 4 + col] += this.data[row * 4 + k] * other.data[k * 4 + col];
        }
      }
    }
    return result;
  }

  public Vec3 transform(Vec3 v) {
    float x = v.x*data[0] + v.y*data[1] + v.z*data[2] + data[3];
    float y = v.x*data[4] + v.y*data[5] + v.z*data[6] + data[7];
    float z = v.x*data[8] + v.y*data[9] + v.z*data[10] + data[11];
    float w = v.x*data[12] + v.y*data[13] + v.z*data[14] + data[15];

    // Convert back from homogeneous coordinates to 3D coordinates.
    if (w != 1.0 && w != 0.0) {
      x /= w;
      y /= w;
      z /= w;
    }

    return new Vec3(x, y, z);
  }

  // Transform the direction without translation.
  public Vec3 transformDirection(Vec3 v) {
    return new Vec3(
      v.x*data[0] + v.y*data[1] + v.z*data[2],
      v.x*data[4] + v.y*data[5] + v.z*data[6],
      v.x*data[8] + v.y*data[9] + v.z*data[10]
    );
  }

  public Ray transformRay(Ray ray) { return new Ray(transform(ray.origin), transformDirection(ray.direction), ray.time); }

  public static Mat4 translate(float tx, float ty, float tz) {
    Mat4 result = new Mat4();
    result.data[3] = tx;
    result.data[7] = ty;
    result.data[11] = tz;
    return result;
  }
  public static Mat4 translate(Vec3 v) { return translate(v.x, v.y, v.z); }

  public static Mat4 scale(float sx, float sy, float sz) {
    Mat4 result = new Mat4();
    result.data[0] = sx;
    result.data[5] = sy;
    result.data[10] = sz;
    return result;
  }

  public static Mat4 rotateX(float angle) {
    final float cos = (float)Math.cos(angle), sin = (float)Math.sin(angle);
    Mat4 result = new Mat4();
    result.data[5] = cos;
    result.data[6] = -sin;
    result.data[9] = sin;
    result.data[10] = cos;
    return result;
  }

  public static Mat4 rotateY(float angle) {
    final float cos = (float)Math.cos(angle), sin = (float)Math.sin(angle);
    Mat4 result = new Mat4();
    result.data[0] = cos;
    result.data[2] = sin;
    result.data[8] = -sin;
    result.data[10] = cos;
    return result;
  }

  public static Mat4 rotateZ(float angle) {
    final float cos = (float)Math.cos(angle), sin = (float)Math.sin(angle);
    Mat4 result = new Mat4();
    result.data[0] = cos;
    result.data[1] = -sin;
    result.data[4] = sin;
    result.data[5] = cos;
    return result;
  }

  // (For now), assumes only one axis is enabled.
  public static Mat4 rotate(float angle, boolean x, boolean y, boolean z) {
    if (x) return rotateX(angle);
    if (y) return rotateY(angle);
    return rotateZ(angle);
  }
}
