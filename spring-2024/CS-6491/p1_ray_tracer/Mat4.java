/**
4x4 matrix
Instances are immutable.
Rotation angles are in radians.
*/
public class Mat4 {
  private final float[] data;

  public Mat4() {
    data = new float[16];
    // Initialize identity matrix.
    for (int i = 0; i < 16; i += 5) data[i] = 1.f;
  }

  // Copy constructor
  public Mat4(Mat4 other) { this.data = other.data.clone(); }

  public Mat4 multiply(Mat4 other) {
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

  public Vec3 apply(Vec3 v) {
    float x = v.x * data[0] + v.y * data[1] + v.z * data[2] + data[3];
    float y = v.x * data[4] + v.y * data[5] + v.z * data[6] + data[7];
    float z = v.x * data[8] + v.y * data[9] + v.z * data[10] + data[11];
    float w = v.x * data[12] + v.y * data[13] + v.z * data[14] + data[15];

    // Convert back from homogeneous coordinates to 3D coordinates
    if (w != 1.0 && w != 0.0) {
      x /= w;
      y /= w;
      z /= w;
    }

    return new Vec3(x, y, z);
  }

  public static Mat4 translate(float tx, float ty, float tz) {
    Mat4 result = new Mat4();
    result.data[3] = tx;
    result.data[7] = ty;
    result.data[11] = tz;
    return result;
  }

  public static Mat4 scale(float sx, float sy, float sz) {
    Mat4 result = new Mat4();
    result.data[0] = sx;
    result.data[5] = sy;
    result.data[10] = sz;
    return result;
  }

  public static Mat4 rotateX(float angle) {
    Mat4 result = new Mat4();
    final float cos = (float)Math.cos(angle), sin = (float)Math.sin(angle);
    result.data[5] = cos;
    result.data[6] = -sin;
    result.data[9] = sin;
    result.data[10] = cos;
    return result;
  }

  public static Mat4 rotateY(float angle) {
    Mat4 result = new Mat4();
    final float cos = (float)Math.cos(angle), sin = (float)Math.sin(angle);
    result.data[0] = cos;
    result.data[2] = sin;
    result.data[8] = -sin;
    result.data[10] = cos;
    return result;
  }

  public static Mat4 rotateZ(float angle) {
    Mat4 result = new Mat4();
    final float cos = (float)Math.cos(angle), sin = (float)Math.sin(angle);
    result.data[0] = cos;
    result.data[1] = -sin;
    result.data[4] = sin;
    result.data[5] = cos;
    return result;
  }

  // Applies the rotation angle to the provided enabled axis _in order x, y, z_.
  public static Mat4 rotate(float angle, boolean x, boolean y, boolean z) {
    Mat4 result = new Mat4();
    if (x) result = result.rotateX(angle);
    if (y) result = result.rotateY(angle);
    if (z) result = result.rotateZ(angle);
    return result;
  }
}
