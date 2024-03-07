class Vec2 {
  final float x, y;

  Vec2() { this.x = this.y = 0; }
  Vec2(float x, float y) {
    this.x = x;
    this.y = y;
  }

  Vec2 add(Vec2 o) { return new Vec2(x + o.x, y + o.y); }
  Vec2 sub(Vec2 o) { return new Vec2(x - o.x, y - o.y); }
  Vec2 mult(float scalar) { return new Vec2(x*scalar, y*scalar); }
  Vec2 div(float scalar) { return new Vec2(x/scalar, y/scalar); } 
  Vec2 flipX() { return new Vec2(-this.x, this.y); }
  Vec2 flipY() { return new Vec2(this.x, -this.y); }

  float length() { return (float)Math.sqrt(x*x + y*y); }
}
