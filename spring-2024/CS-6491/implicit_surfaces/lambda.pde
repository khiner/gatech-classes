// Lambda expressions for implicit functions
//
// See the "a_sphere" lambda expression for an example of defining an implicit function

import java.lang.FunctionalInterface;

// this is a functional interface that will let us define an implicit function
@FunctionalInterface
interface ImplicitInterface {
  // abstract method that takes (x, y, z) and returns a float
  float getValue(float x, float y, float z);
}

// Implicit function for a sphere at the origin.
//
// This may look like a function definition, but it is a lambda expression that we are
// storing in the variable "a_sphere" using =. Note the -> and the semi-colon after the last }

ImplicitInterface a_sphere = (x, y, z) -> {
  float d = sqrt(x*x + y*y + z*z);
  return (d);
};
