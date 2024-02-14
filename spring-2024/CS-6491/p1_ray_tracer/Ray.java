class Ray {
  final Vec3 origin, direction;
  
  Ray(Vec3 origin, Vec3 direction) {
    this.origin = origin;
    this.direction = direction;
  }

  // Interpolate along the ray to find the point at distance `t`.
  Vec3 interp(float t) { return origin.add(direction.mult(t)); }
}
