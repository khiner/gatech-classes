class PointLight extends Light {
  final Vec3 position;

  PointLight(Vec3 position, Color c) {
    super(c);
    this.position = position;
  }

  Vec3 samplePoint() { return position; }
}
