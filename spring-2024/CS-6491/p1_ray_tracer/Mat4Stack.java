import java.util.Stack;

/*
The stack is initialized with the identity matrix and is always non-empty.
Rotation angles are in radians.
*/
public class Mat4Stack {
  private final Stack<Mat4> stack;

  public Mat4Stack() {
    stack = new Stack<>();
    stack.push(new Mat4());
  }

  public Mat4 top() { return stack.peek(); }

  public void push() {
    stack.push(new Mat4(top()));
  }
  public void pop() {
    if (stack.size() > 1) stack.pop();
  }
  
  // Apply the provided transform to the top of the stack.
  public void apply(Mat4 transform) {
    stack.push(stack.pop().mult(transform));
  }
}
