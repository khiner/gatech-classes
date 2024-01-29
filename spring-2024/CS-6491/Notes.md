# Graphics notes

## Fri Jan 12

### Eye rays

x(t) = x_0 + t(x_1 - x_0)

Translate x \in [0,w], y \in [0,h] -> x', y' \in [-k,k] (the view plane)

x' = (x - w/2)*(2k/w)
y' = (y - h/2)*(-2k/h)

### Ray triangle intersection

Implicit plane equation

f(x,y,z) = ax + by + cz + d ?= 0

Given triangle: P,Q,R (3d points)

PQ = Q-P, PR = R-P
surface normal:
N = \frac{PQ \cross PR}{\|PQ \cross PR\|}
N is a unit length vector, perp. to plane of P,Q,R

N = (a,b,c). But what aboud d?
Substitute P = (x,y,z) into implicit plane equation f and solve for d.

t(a dx + b dy + c dz) + ax_0 + b_y0 + cz_0 + d = 0
t = \frac{-(ax_0 + by_0 + cz_0 + d)}{a dx + b dy + c dz}

If num == den, ray parallel w/ tri and it's a miss.

### Point-in-triangle

Triangle: P,Q,R

Half-plane test (2D):
Extend the triangle's edges infinitely.
Assign +/- to given point based on which side of each line it's on, and if all 3 are +, it's in triangle.

2d - which side of line thru AB is P?
V = AP \cross AB = (0,0,V_z)
if V_z == 0, P is on the line AB. If V_z > 0, P is on one, else on the other

side(A,B,P) -> (AP \cross AB)\_z
point_in_triangle(P,Q,R,p) -> side(P,Q,p) == side(Q,R,p) == side(R,P,p)

3D:
V _ N is + or - (where N is surface normal of triangle)
Go around counterclockwise (or clockwise)
side(A,B,P,N) = sign of N _ (AP \cross AB)
side(B,C,P,N) = ...
side(C,A,P,N) = ...
All must be equal sign.

## Fri Jan 19

Scene graphs:
    Scene 1: (TS) -> P
    P'' = T * S * P
    Scene 2: ST -> P
    P'' = S * T * P

Box with extent 1, in 3D centered at origin
w.t. move it by T = translate(0,0,-10)

push()
translate(0,0,-10)
box()
pop()

matrix stack
top of stack = current transformation matrix C
initially, C = I

## Mon Jan 22

Bounding volumes - computational savings for ray misses
```
class Box:
    expandBy(Point)
    expandBy(Box)
    expandByEpsilon() // expand to next machine-precision floats

class Polygon:
    List<Triangle>

    Box getBbox()
```

Intersect ray with axis-aligned plane

Define plane p as all solutions for x = x_p (only need a single number to describe an implicit equation for an axis-aligned plane)
ray: a + tb
x(t) = a_x + tb_x
(solve for x_p):
t = (x_p - a_x)/b_x // point of intersection with plane along ray

For axis-aligned bbox intersection, just do this multiple times.

in 2d, find x/y-intervals
[t_min_x_intersect, t_max_x_intersect]
[t_min_y_intersect, t_max_y_intersect]
If intervals overlap (or if exactly one interval is empty), ray _does_ intersect plane.

Bounding volume hierarchies
```
class Node {
    Box bbox;
    List<Node> children;
    List<Object> objects;
}

traverse(Ray r, Node n) {
    if (r misses n's bbox) return null;
    if (n is leaf) return (intersect r with n's objects);
    for (Node child c : n.children) {
        hit = traverse(r, c);
        if (hit) return hit;
    }
    return null;
}
```

Build hierarchy top-down:
- start: bbox enclosing all objects
- pick a point along an axis-aligned plane and split the list into 2 bboxs (child nodes)
- stop when objecct list is "short"

## Wed Jan 24

centroid of bbox: (cx,cy,cz) - center of mass
cx = (xmin + xmax)/2

choose split axis:
find min + max range of centroids
pick split axis with longest range

how to split the axis:
- option 1: split on middle
- option 2: split based on even counts:
    count to (n/2)th bounding box

3D grid: uniform spatial subdivision
Each cell has a list of objects inside it
The ray intersects the 3D grid at regular distances _in any single axis_.

When we intersect a shape that's at least partially in a given cell, we need to see if the intersection point is actually in the cell.

We can use the "axis-aligned bbox intersection" above to check.

Intersection "y" points, t_i, at each cell x axis x_i:
t_i = t_(i-1) + h/b_x,
where h is the horizontal grid cell size

3d grid traversal:
step separately along each axis

Per-axis values (x,y,z)
___
delta_t (float): change in t that goes to next cell
step (int): 1 or -1, amount to add to go to next cell
out (int): cell value indicating we've left the grid
pos (int): current cell we are checking
next_crossing_t (float): what value of t means we have crossed to next cell

Step to next cell (move by one cell)
___
while ...
    Determine smallest next_crossing_t, call this step_axis (0,1,2)
    pos[step_axis] += step[step_axis]
    if (pos[step_axis] == out[step_axis]) break;
    next_crossing_t[step_axis] += delta_t[step_axis]

(we call it t because of ray: a + tb)

## Fri Jan 26

Reflection, refraction, shadows

ray(t) = P + tD

P (shadow ray)
visible(P,Q) == visible(Q,P)

Reflection:
    - Angle of incoming ray is reflected about the surface normal of the reflecting surface
    - Effectively, we move the origin of our incoming ray (eye) to the intersection point of the reflecting surface (virtual eye) to create the reflected ray

Create reflected ray:
    Assume `N` (normal of reflective surface) is unit-length.
    Vector `v` _pointing at_ the eye.
    Project `v` onto `N`: `N(N \dot v)`
    Multiply by 2 (double the length), then subtrace `v`
    Final reflected direction: `2N(N \dot v) - v`

When to top bouncing off of reflective surfaces in the scene?
    surface color: `c = diffuse + k_{refl}c_{refl}`
    `k`: how perfect a mirror, `k <= 1`,
    `refl`: label for reflective surface
    Recursive! So `k` gets smaller and smaller...
    End recursion
        (1) when contribution of ray `c` is small (e.g. `< 0.01`)
        (2) simpler way: stop when ray depth is too high (e.g. 10-20 reflections)

Refraction = Transmission:
Transmitted rays: (potentially bent) rays within an object with different material properties than outside the object
    Main cause: Speed of light is different in different materials
    Light travels fastest through a vacuum
    Index of refraction of material ($n$ or $\eta$) = (SOL through vacuum) / (SOL through material)
        Always >= 1 by definition
    Snell's law: $\frac{\sin{\theta_i}}{\sin{\theta_t}}$
    $\theta_i$ incident angle, $\theta_t$ transmitted angle
    Beyond a "critical angle", all light is reflected at the surface of the object _back into_ the object
    Used in e.g. fiber optics, to keep light inside cable

Fresnel effect:
    glass, plastic, ceramic, paper, ... "dialectrics"
    `n = 1` (air)
    `n_t >> 1`
    Instead of constant $k$,
    $$k = R(\theta) = R_0 + (1-R_0)(1-\cos(\theta))^5$$
    $$R_0 = (\frac{n_t - 1}{n_t + 1})$$, usually < 1
    This is the Schlick approximation to Fresnel
    steep angle -> lower refl, shallow angle -> higher refl.
    $$k_{trans} = T(\theta) = 1 - R(\theta)$$

## Mon Jan 29

### Surface Generation

Implicit surfaces
An implicit surface is the set of zeros of a function of three (in 3D) variables.

Implicit function: tell whether point is in, out, or on a surface
Often use distances or dist^2

Ex. f(p) = r^2 - \|p-c\| p = (x,y,z) c = (x_c,y_c,z_c)
(we can flip the sign here based on convention, whether we consider pos. values inside (as above) vs. outside the surface)

circle or sphere (2d or 3d) = r^2 - [(x-x_c)^2 + (y-y_c)^2 + (z-z_c)^2]

Basis functions / falloff functions:
g(d) = e^{-d^2} Gaussian function (used by James Blinn, so sometimes called Blinn function)

d(d) = (1-d^2)^3 (Wyvill - this actually crosses zero, unlike Blinn. We only use the pos. part)

Use with distance functions (to points):
    f_1(p) = g(\|p-c_1\|)
    f_2(p) = g(\|p - c_2\|)

Form sum:
    f(p) = f_1(p) + f_2(p) - T,
    where T is a scalar

Distance to line segment

D = P2 - P1
Point Q
Project Q onto D: V = Q-P1, t = V \dot D / |D|^2
Closest point P on segment D: {t < 0 -> P1, t > 1 -> P2, else t in [0,1] -> P1 + tD}

d(Q) = \|P-Q\|
blobby form: f(Q) = g(d(Q))
