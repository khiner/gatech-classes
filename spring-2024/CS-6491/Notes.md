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
- stop when object list is "short"

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

## Wed Jan 31

### Ray tracing instancing

r(t) = a + tb

Transform polygon by C : transform n vertices (3 for triangle), _then_ ray trace

How to ray trace object with non-uniform scaling?
    Transform the _ray_ insterad of the object! "A ray transformed is just another ray."

Canonical sphere at the origin
Non uniform transformation C
Ray trace with ray r - transform the _ray_ by C^{-1} (a, b) -> (C^{-1}a, C^{-1}b)
(Transformed space C to "canonical" space C^{-1})

C^{-1}a = C^{-1}[b_x b_y b_z 0] (col vec)

transform the surface normal n (in canonical space)
    n -> (C^{-1}^T)n

So: DON'T transform every vertex!

### Iso-surface extraction

Render implicits -> (ray trace | convert to polygons & render polys)

2D marching squares
Sample grid. Estimate zeros by linearly interpolating between points that switch sign.

3D: Marching cubes: 8 corneers of cube, 12 edges

### Dual contouring

Estimate edge inersection points, and the gradient of f at those points,
and find interior points at the intersection of the gradient lines (instead
of putting all points on the grid lines).

Better at capturing curves.

## Fri Feb 2

Implicit surfaces manipulation: good section in textbook

Distance function: d(p)
Blob function: g(d) = {
    (1-d^2)^3    d <= 1
    0    d > 1
}

twist(x,y,z) = {
    x
    y * cos(x) - z * sin(x)
    y * sin(x) + z * cos(x)
} (use x like an angle)

Taper

scale shape by k -> divide(x,y,z) by k
...

t(x) = {
    0    x < x_min>
    (x-x_min)/(x_max-x_min)    x_min <= x <= x_max
    1    x > x_max
}

k(x) = (1-t(x))k_1 + t(x)k_2

taper(x,y,z) = {
    x
    y / k(x)
    z / k(x)
}

Common warps: twist,taper,bend - due to Alan Barr siggraph 1984

bend: see fundamentals ch 22


### Constructive Solid Geometry (CSG)

boolean combinations of shapes

intersect(a,b) (and)
union(a,b) (or)

use `min` for implicit intersection (minimum of two functions) and `max` for union

intersect(A(p), B(p)) = min(A(p),B(p)) => A and B
union(A(p), B(p)) = max(A(p),B(p)) => A or B
subtract(A(p), B(p)) = min(A(p), 2 * T - B(p)) => A and not B

### 3d morphing

Linear interp between functions:
morph(A,B,t) = (1-t)A(p) + tB(p)


## Mon Feb 5

### Ray tracing implicit surfaces

Simple shapes can be ray traced analytically (plane, sphere, ellipsoid, cylinder, cone)

p(t) = a + tb: assume unit length b, t = distance along the ray

Ray/implicit inersection = root finding

f(p(t)) = f(a + tb): 1D function in t

Main approach: start at p = ray origin = p(0)
take small steps along ray until we find f(p(t)) = 0

Strategy 1 (Ray Marching): Step by constant delta-t along ray (slow)
    small t -> more accurate intersections
    large t -> faster, looks worse
Then use binary search to find the zero after we find transition from neg->pos

Strategy 2 (Sphere Tracing): Step based on knowing distance to surface (faster)
Want: take large steps if we are far from surface, smal when near surface
How do you know dist. from p(t) to surface?
    1. you have a sampled distance field (3d grid of distance values)
    2. when gradient magnitude (of scalar field f(p)) is bounded
        ($\|\Delta f(p) \| < k$ - implied by Lipschitz condition$)
        Suppose gradient magnitude (slope) is at most 1 (1-Lipschitz)
        Then, can take step $\Delta t = \|f(p)\|$ (step based on _function value_).


## Wed. Feb 7

Mesh: surface made of poygons
Polyhedron has polygon foces. Made of vertices and edges.

Manifold: all mesh regions look "like a plane"
    (not all meshes are manifold)

Manifold vertex: Can circumnavigate around the vertex, hitting _every_ surrounding face ("fan" of faces)
Non-manifold vertex: "Fall off the edge" when trying to walk around it

Manifold edge: Has exactly two neighboring faces
Non-manifold edge: Has more than two neighboring faces
Boundary edge: Only one polygon neighbor face (not manifold edge, but not so bad either)

### Common mesh operations
- Laplacian smoothing
new_v = (v1+v2+v3+v4+v5)/5
    - for every v
Smooths surface, may shrink it, improved triangle shape
- Face subdivision
Add more faces in the middle of existing faces
Used to make geodesic spheres
- Triangulation
Make all faces triangles
Why tris? Always planar! Fixed size data structure. Interpolation is well defined (barycentric interpolation).

### Mesh representations

Polygon soup
    begin vertex vertex vertex end
    ...
- Creates large files due to repeated vertices
    - used in .stl files (used for 3d printing)

Indexed mesh
- Common file format (.obj, .ply)
- List of vertices, and list of faces w/ vertex indices


## Fri Feb 9

### Mesh connectivity/adjacency/neighbors

More mesh ops: simplify, hole filling, ...

4 main kinds of adjacency:
* face->vertices (store in face)
* face->adjacent faces
* vertex->adjacent faces
* vertex->vertices (connected by an edge)

Like indexed mesh, store vertices of face in CCW order

Other adjacency from: directed edge: Bridge between faces + vertices

Directed edge: edge of face, from v1 to v2 (order matters!)
e.g. two adjacent tris, f1 & f2. v1 on bottom, v2 on top. f1 left tri, f2 right tri.
    e1 = (f1,v1,v2), e2 = (f2,v2,v1)

Find + store opposite edge info: map
    key: v1, v2
    value: face, vertex index of v1 <- edge info

1) store all (face,vertex index) mesh for all faces -> map
2) find opposite (v2,v1) for each edge at face, store this in the face

class Vertex {
    float x,y,z;
    Face one_face;
}

class Edge {
    Face f;
    int v_index;
}

class Face {
    Vertex verts[];
    Edge opposite_edges[];
}

mesh:
    list of vertices
    list of faces

Edge ops:
    next(e) - next edge ccw in face
    prev(e) - previous edge ccw in face
    opp(e) - opposite edge (stored in face)
    swing(e) - next edge cw around start vertex (next(opp(e)))


vertex->face/vertex->vertices: pseudo-code
    e_start = edge_from_vertex(Vertex v, Face v.one_face);
    e = e_start;
    do {
        // do something with face of e
        e = swing(e);
    } while(e != e_start);

similar to half-edge


## Mon Feb 12

Regular planar tilings: Fill the plane with regular polygons
all one kind, no overlaps, no gaps

Square, triangular, and hexagonal tiling are all regular planar tilings.

valence of vertex/face: num faces surrounding it

Dual: swaps faces and vertices
    put a vertex in center of each face.
    rotate each edge separating original faces 90 degrees to connect the center vertices

Square tiling is _self-dual_.

Dual of triangular tiling is hexagonal.

Dual switches vertex valence with num face sides

Platonic solids
- Tetrahedron: face -> 3 sides, vertex -> valence 3
    - Tetrahedron is _self-dual_.
- Octahedron: face -> 3 sides, vertex -> valence 4
    - Dual of is a cube
- Cube: face -> 6 sides, vertex -> valence 8
    - Dual is Octahedron
- Icosahedron: face -> 3 sides, vertex -> valence 5
    - Dual is dodecahedron
- Dodecahedron: face -> 5 sides, vertex -> valence 3? (double check)

Name|Vertices|Faces|Edges
--|--|--|--
Tetrahedron|4|4|6
Octahedron|6|8|12
Cube|8|6|12
Icosahedron|12|20|30
Dodecahedron|20|12|30

V+F=E+2 (for manifold, connected, genus zero (no holes))

_Not_ platonic: Prism
Start out with n-gon. Place a copy below, and connect them.

Anti-prism is a prism with one side rotated.


## Wed Feb 14

Mesh properties

Property|# dimensions|count symbol
--|--|--
Vertices|0|V
Edges|1|E
Faces/Trianges|2|F/T
Cells/Volumes|3|_

A mesh is _oriented_ if all edges have opposite order (all faces in same ccw/cw ordering)

Counterexample: Cannot orient a Klein bottle!

A mesh is _connected_ if there is a path between any 2 vertices.
A manifold mesh has all edges and vertices manifold
A shell is one connected component of a mesh

A handle is a hold through a mesh (sphere has no handles, torus has one handle)
Num handles is _genus_ of the surface.

Let G be genus of surface (num holes)
Let S be num shells of surface
General Euler-Poincare equation: V - E + F = 2(S - G)
Special cases:
    When S = 1: V - E + F = 2 - 2G = $\Chi$ (Euler characteristic of a surface)
    Triangle mesh (manifold w/o boundary, G = 0): V + T = E + 2.
    Also, **2E = 3T**, **2V = T + 4**, and **3V = E + 6**

genus:0 1 2 3 $\to$ $\Chi$: 2 0 -2 -4

- V: X
- T: X X
- E: X X X

Let $\bar{D}$ be average degree (valence) of vertices.

Create $\bar{D}$ edges for V vertices, creates double the edges:
    $\bar{D}V = 2E \implies \bar{D} = 2(3V - 6)/V \approx 6$, for large V
So, for triangle meshes, average valence $\approx$ 6

Quadrilaterals: $\bar{D} = 4$


