class Color {
  final float r, g, b; // In range [0, 1]

  Color(float r, float g, float b) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  
  Color add(Color other) { return new Color(r + other.r, g + other.g, b + other.b); }
  Color add(float scalar) { return new Color(r + scalar, g + scalar, b + scalar); }

  Color mult(Color other) { return new Color(r * other.r, g * other.g, b * other.b); }
  Color mult(float scalar) { return new Color(r * scalar, g * scalar, b * scalar); }

  Color div(Color other) { return new Color(r / other.r, g / other.g, b / other.b); }
  Color div(float scalar) { return new Color(r / scalar, g / scalar, b / scalar); }
}
