// Axis-aligned bounding box
class BBox extends Geometry {
  final Vec3 min, max;

  BBox(Vec3 min, Vec3 max) {
    this.min = min;
    this.max = max;
  }

  static final float MAX = Float.MAX_VALUE;
  static final BBox EMPTY = new BBox(new Vec3(MAX,MAX,MAX), new Vec3(-MAX,-MAX,-MAX));

  int maxAxis() { return max.sub(min).maxAxis(); }
  Vec3 center() { return min.add(max).div(2); }
  BBox union(BBox o) { return new BBox(Vec3.min(min, o.min), Vec3.max(max, o.max)); }

  // Normal based on the closest axis-aligned plane.
  Vec3 normal(Vec3 p) {
    final float eps = 1e-4f;

    if (Math.abs(p.x - min.x) < eps) return new Vec3(-1, 0, 0);
    if (Math.abs(p.x - max.x) < eps) return new Vec3(1, 0, 0);
    if (Math.abs(p.y - min.y) < eps) return new Vec3(0, -1, 0);
    if (Math.abs(p.y - max.y) < eps) return new Vec3(0, 1, 0);
    if (Math.abs(p.z - min.z) < eps) return new Vec3(0, 0, 1);
    if (Math.abs(p.z - max.z) < eps) return new Vec3(0, 0, -1);

    return null;
  }

  /*
  * Based on PBRTv4's
  * [Bounds3::IntersectP](https://github.com/mmp/pbrt-v4/blob/39e01e61f8de07b99859df04b271a02a53d9aeb2/src/pbrt/util/vecmath.h#L1546C1-L1572C1)
  */
  Float intersectT(Ray ray) {
    final float eps = 1e-5f;

    float t0 = 0, t1 = MAX;
    for (int i = 0; i < 3; ++i) {
      final float d = ray.direction.at(i), o = ray.origin.at(i);
      final float min_val = min.at(i), max_val = max.at(i);
      // This will be Infinity when d == 0, and subsequent calculations work out correctly.
      // Even though ops involving Inf/NaN can be slow, in my testing this is faster than introducing branching.
      final float invD = 1.f / d;

      float tNear = (min_val - o) * invD, tFar = (max_val - o) * invD;
      if (tNear > tFar) {
        float temp = tNear;
        tNear = tFar;
        tFar = temp;
      }

      tFar *= 1 + 2*eps;
      t0 = Math.max(t0, tNear);
      t1 = Math.min(t1, tFar);
      if (t0 > t1) return null;
    }

    return t0;
  }

  Intersection intersect(Ray ray) {
    final Float t = intersectT(ray);
    if (t == null || t == 0) return null;


    final Vec3 norm = normal(ray.interp(t));
    if (norm == null) throw new RuntimeException("Ray intersected BBox but could not determine the normal based on the intersection point.");

    return new Intersection(ray, t, norm);
  }

  boolean intersects(Ray ray) { return intersectT(ray) != null; }

  BBox getBBox() { return this; }
}
