public class Vec3 {
  float x, y, z;

  public Vec3(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public Vec3 add(Vec3 o) { return new Vec3(this.x + o.x, this.y + o.y, this.z + o.z); }
  public Vec3 add(float v) { return new Vec3(this.x + v, this.y + v, this.z + v); }

  public Vec3 sub(Vec3 o) { return new Vec3(this.x - o.x, this.y - o.y, this.z - o.z); }
  public Vec3 sub(float v) { return new Vec3(this.x - v, this.y - v, this.z - v); }

  public Vec3 mult(Vec3 o) { return new Vec3(this.x * o.x, this.y * o.y, this.z * o.z); }
  public Vec3 mult(float v) { return new Vec3(this.x * v, this.y * v, this.z * v); }

  // In-place
  public void flip() {
    this.x *= -1;
    this.y *= -1;
    this.z *= -1;
  }
    
  public float dot(Vec3 o) { return this.x * o.x + this.y * o.y + this.z * o.z; }

  public Vec3 cross(Vec3 o) {
    return new Vec3(
      this.y * o.z - this.z * o.y,
      this.z * o.x - this.x * o.z,
      this.x * o.y - this.y * o.x
    );
  }

  public Vec3 normalize() {
    final float mag = length();
    return mag == 0 ? new Vec3(0, 0, 0) : new Vec3(x / mag, y / mag, z / mag);
  }

  public float length() { return (float)Math.sqrt(x*x + y*y + z*z); }
}
