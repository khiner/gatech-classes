public class Vec3 {
  final float x, y, z;

  public Vec3(float x, float y, float z) {
    this.x = x;
    this.y = y;
    this.z = z;
  }

  public float at(int i) { return i == 0 ? x : i == 1 ? y : z; }
  public int maxAxis() {
    final Vec3 absv = abs();
    if (absv.x >= absv.y && absv.x >= absv.z) return 0;
    if (absv.y >= absv.x && absv.y >= absv.z) return 1;
    return 2;
  }

  public Vec3 add(Vec3 o) { return new Vec3(this.x + o.x, this.y + o.y, this.z + o.z); }
  public Vec3 add(float v) { return new Vec3(this.x + v, this.y + v, this.z + v); }

  public Vec3 sub(Vec3 o) { return new Vec3(this.x - o.x, this.y - o.y, this.z - o.z); }
  public Vec3 sub(float v) { return new Vec3(this.x - v, this.y - v, this.z - v); }

  public Vec3 mult(Vec3 o) { return new Vec3(this.x * o.x, this.y * o.y, this.z * o.z); }
  public Vec3 mult(float v) { return new Vec3(this.x * v, this.y * v, this.z * v); }

  public Vec3 div(Vec3 o) { return new Vec3(this.x / o.x, this.y / o.y, this.z / o.z); }
  public Vec3 div(float v) { return v != 0 ? new Vec3(this.x / v, this.y / v, this.z / v) : new Vec3(0, 0, 0); }

  public Vec3 flip() { return new Vec3(-this.x, -this.y, -this.z); }
  public Vec3 abs() { return new Vec3(Math.abs(x), Math.abs(y), Math.abs(z)); }
    
  public float dot(Vec3 o) { return this.x * o.x + this.y * o.y + this.z * o.z; }

  public Vec3 cross(Vec3 o) {
    return new Vec3(
      this.y*o.z - this.z*o.y,
      this.z*o.x - this.x*o.z,
      this.x*o.y - this.y*o.x
    );
  }

  public Vec3 normalize() { return div(length()); }

  public float length() { return (float)Math.sqrt(x*x + y*y + z*z); }

  public static Vec3 min(Vec3 a, Vec3 b) { return new Vec3(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.min(a.z, b.z)); }
  public static Vec3 max(Vec3 a, Vec3 b) { return new Vec3(Math.max(a.x, b.x), Math.max(a.y, b.y), Math.max(a.z, b.z)); }
}
