import java.util.Stack;

/*
The stack is initialized with the identity matrix is always non-empty.
Rotation angles are in radians.
*/
public class Mat4Stack {
  private final Stack<Mat4> stack;

  public Mat4Stack() {
    stack = new Stack<>();
    stack.push(new Mat4()); // Initialize with identity matrix.
  }

  public void push() {
    stack.push(new Mat4(stack.peek()));
  }
  public void pop() {
    if (stack.size() > 1) stack.pop();
  }
  
  // Apply the provided transform to the top of the stack and replace it.
  public void apply(Mat4 transform) {
    stack.push(transform.multiply(stack.pop()));
  }
  
  public Mat4 top() { return stack.peek(); }

  // Convenience methods to apply specific transformations to the top of the stack.
  public void translate(float tx, float ty, float tz) { apply(Mat4.translate(tx, ty, tz)); }
  public void scale(float sx, float sy, float sz) { apply(Mat4.scale(sx, sy, sz)); }
  public void rotateX(float angle) { apply(Mat4.rotateX(angle)); }
  public void rotateY(float angle) { apply(Mat4.rotateY(angle)); }
  public void rotateZ(float angle) { apply(Mat4.rotateZ(angle)); }
}
