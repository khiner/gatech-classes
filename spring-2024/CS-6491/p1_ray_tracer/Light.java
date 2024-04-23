abstract class Light {
  final Color c;

  Light(Color c) {
    this.c = c;
  }

  abstract Vec3 samplePoint();
}
