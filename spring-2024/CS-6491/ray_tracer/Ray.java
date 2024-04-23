class Ray {
  final Vec3 origin, direction;
  final float time; // Time within the frame -> [0,1]
  
  Ray(Vec3 origin, Vec3 direction, float time) {
    this.origin = origin;
    this.direction = direction;
    this.time = time;
  }

  // Interpolate along the ray to find the point at distance `t`.
  Vec3 interp(float t) { return origin.add(direction.mult(t)); }
}
