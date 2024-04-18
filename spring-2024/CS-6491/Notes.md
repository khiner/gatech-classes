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
    0    x < x_min
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
    - put a vertex in center of each face.
    - rotate each edge separating original faces 90 degrees to connect the center vertices

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
- F: X X
- E: X X X

Let $\bar{D}$ be average degree (valence) of vertices.

Create $\bar{D}$ edges for V vertices, creates double the edges:
    $\bar{D}V = 2E \implies \bar{D} = 2(3V - 6)/V \approx 6$, for large V
So, for triangle meshes, average valence $\approx$ 6

Quadrilaterals: $\bar{D} = 4$


## Mon Feb 19

Line segments

Q(t) = p1 + t(p2-p1) = (1-t)p1 + tp2
here, _weights_ l1(t) = 1-t, l2(t) = t

l1(t)+l2(t)=1, non-negative if 0 <= t <= 1

Q(t) is convex comb. of p1 and p2

### Quadratic (degree 2) Bezier curve

3 control points: p1,p2,p3

t=0.7, subdivide p1,p2 and p2,p3 by 0.7, draw a line between those points, and subdivide _that_ line at 0.7 - that point is Q(t=0.7)

- Q(t) = S1(t)P1 + S2(t)P2 + S3(t)P3
- S1(t) = (1-t)^2
- S2(t) = 2t(1-t)
- S3(t) = t^2
- $\implies$ S1+S2+S3=1

### Cubic (degree 3) Bezier curve

4 control points: p1,p2,p3,p4

1) independent tangent control at endpoints
2) curve inside convex hull of control points
3) curve interpolates the points p1 and p4
4) curve is tangent to p1-p2 at t=0, and tangent to p3-p4 at t=1

subdivide each line by t, then you have two lines then you do the same subdivision as Quadratic Bezier to get Q(t) on a line between those two lines

- Q(t) = B1(t)P1 + B2(t)P2 + B3(t)P3 + B4(t)P4
- B1(t) = (1-t)^3
- B2(t) = 3t(1-t)^2
- B3(t) = 3t^2(1-t)
- B4(t) = t^3
- $\implies$ B1+B2+B3+B4=1


## Wed Feb 21

### B-splines

- curve does not pass through control points
- cubic b-splines, c2 continuous
- sliding windows of points (4 at a time)

- S1(t) = 1/6 (1-t)^3
- S2(t) = 1/6 (3t^3 - 6t^2 + 4)
- S3(t) = 1/6 (-3t^2 + 3t^2 + 3t + 1)
- S4(t) = 1/6 t^3
- $\implies$ S1+S2+S3+S4=1

### Subdivision Surfaces

Given: base mesh (polygons)
Result: smooth surface

First published: 1978 Catmull-Clark, Doo-Sabin

#### Loop Subdivision (Triangles)

(not loops! invented by Charles Loop, lol)
C2 continuous nearly everywhere

Subdivision:
1) compute locs of new verts (associated with edges)
2) move positions of old verts, too
3) make smaller triangles 1->4

New edge assoc. vertex:
v = 3/8 (v1 + v2) + 1/8 (v3 + v4)

k is valence of v

new vertex: v':
v' = (1-k\beta)v + \beta(v1+v2+...+vk)
\beta = 3/(8k) for k > 3, \beta=3/16 if k = 3

(1-k/beta) + /beta k = 1

only c1 continuous (at the limit) at "extraordinary points" (non-6 valence)


## Fri Feb 23

Butterfly interpolation

...

### Catmul-Clark

New face vertices for polygon faces

quad: v = 1/4(v1+v2+v3+v4)
k-gon: v = (centroid) = \frac{1}{k} \sum_{i=1}^k{v_i}

Let E = 1/k (e1+e2+...+ek)
    F = 1/k (f1+f2+...+fk)

v' = 2E/k + F/k + v(k-3)/k
k is valence of vertex v

Most vertices valence 4

Crease Marking - can "tag" an edge to be sharp

Semi-Sharp Edges: Use sharp edge rule S times,
with S = 0 for no sharpness, S = 1 a little sharper, etc.


## Mon Feb 26

### Smoothing

#### 1D function smoothing

use 2nd derivative to find min+max and use these: decrease/increase function

f_new = f_old + \lambda f'', 0 <= \lambda <= 1

#### 2D curve smoothing

push convex parts in, pull convcave parts out

Approximate curve as points p_i

1st derivative
- t_{i-1} = p_i - p_{i-1}
- t_i = p_{i+1}-p_i

2nd derivative (Laplacian)
\Delta p_i = 1/2 (t_i - t_{i-1}) (defference of tangents)
    = 1/2[(p_{i+1} - p_i) - (p_i - p_{i-1})]
    = 1/2(p_{i+1}+p_{i-1}) - p_i (First term is avg of neighbors)

p_i += \lambda \Delta p_i

#### 3D mesh smoothing

same, w/ \Delta v = v_c - v,
where v_c is center vertex

Don't move vertices one at a time!
(Keep old positions while calculating new ones and replace all at once.)

Fix shrinkage: Alternate shrinking and inflating (Gabriel Taubin)

Repeat:
    v += \lambda \Delta v   \lambda > 0 (shrink)
    v += \mu \Delta v   \mu < 0 (inflate)
    \lambda = 0.6307, \mu = -0.67315

**kesen** - Ke-Sen Huang - siggraph paper archives

## Mon Mar 4

### Radiometry

Reflectance constant of proportionality:
Bidirectional Reflectance Distribution Function (BRDF)
$$\rho(k_i,k_o) = \frac{\text{outgoing light towards}\ k_o}{\text{incoming light from}\ k_i} = \frac{dL(x,k_o)}{L(x,k_i)\cos(\theta)dk_i}$$
Describes how a material reflects light.

Function of 4 variables: $\rho(k_i,k_o), k_i = (\theta_i,\phi_i)$

Hemisphere = $\Omega$
- All possible incoming + outgoing directions in the hemisphere.
- Incoming light at a point $p$ on surface. Hemisphere $\Omega$ centered at $p$.
- Angle $\theta$ around the pole, from the top
- Angle $\phi$ around base

BRDF Laws:
1) Helmholtz reciprocity: $\rho(k_1,k_2) = \rho(k_2,k_1)$
  - photons travel on the same path forward and backward
2) Conservation of energy:
$$\int_{k_i \in \Omega}{\rho(k_i,k_o)cos\theta_i dk_i} \leq 1$$
- Energy out is $\leq$ incoming energy

- Light falling on surface:
    -  $\theta$: Angle of surface relative to incoming light
    - $A$: 1D unit area of surface _perpendicular to light_
    - $D$ = 1D unit area of surface = $A / cos(\theta)$

For many surfaces, $\phi_i$ (incoming angle) can be ignored $\to$ isotropic surfaces
- For isotropic surfaces, $\rho(\theta_i, \theta_o, \phi_o)$ - 3D functions
- Anisotropic surfaces: more complex (brushed metal, velvet, some fabrics, hair)

Easier to show 2D diagrams of BRDFs: (slice of hemisphere)

- Half circle, with zero at top, -90deg left and 90deg right
- Incoming light dir: $\theta_i$
- Reflects off surface in multiple directions ("specular" glossy/fairly shiny surface)
- Diffuse surface: Outgoing light scattered in all directions.
- Perfect mirror surface: All outgoing light reflects in a single, mirrored direction: $\theta_o = -\theta_i$


## Wed Mar 6

### Distributed (Distribution) Ray Tracing

Siggraph 1984: Rob Cook, ...

Soft rendering effects:
    - soft shadows
    - glossy reflection (between perfect mirror and diffuse)
    - motion blur
    - depth of field

Soft shadows, with area light
    - Umbra: Total shadow
    - Penumbra: Partial shadow

Traditional ray tracing:
    - From eye to point light, one path
    - bool visible(P,L)
    - C = (N\dot L)visible(P,L)

Distribution ray tracing:
    - From eye to area light
    - C=1/n sum_i(visible(P,L_i)) Na \dot L_i
    - where to shoot rays from eye?
        - jittered sampling, pick a random point within each grid cell on area light
        - rejection sampling, reject points sampled outside of desired area
        - combine the two.
        - (best: blue noise sampling)

Glossy reflection
    - Don't only bounce as a perfect mirror. Avg. multiple scatters (BRDF)
    - Mirror: C = diffuse + k_refl C_refl
    - Glossy: C = diffuse + avg(k_refl C_R_i)

Motion blur
    - Ray through view plane
    - Perform multiple intersections with a single ray w/ _multiple timestamps_
    - Jitter sample points in time within "shutter" window

DOF:
    - Pinhole camera forms upside-down image from light passing only through the pinhole
    - Let's widen the hole into the box. 2D plane parallel to the lens. Lens bends the rays to a focal point on the film plane. Points _not_ on focal plane don't converge to a focal point. (Circle of Confusion)
    1) where would light hit if it went exactly through the middle of the lens? Call that P.
    2) Pick another point in the lens, not in the middle, and direct a ray to P.
    3) Send many rays from different points on the lens to P, average them all

## Fri Mar 8

### Global illumination

Distribution ray tracing sums over multiple light paths:
- multiple shadow rays
- multiple reflected rays (glossy)

The two above both approximate the _Reflectance Equation_:
$$L(x,k_o) = \int_{k_i \in \Omega}{\rho(k_i,k_o)L(x,k_i) \cos(\theta_i) dk_i},$$
where $L(x,k_o)$ is the outgoing radiance, $k_i \in \Omega, k_i = (\theta_i, \phi_i)$ is all directions around the hemisphere, $\rho(k_i,k_o)$ is the fraction of incoming light that reflects to $k_o$, $L(x,k_i)$ is the incoming light to $x$ from $k_i$, $\cos(\theta_i)$ is the area-weighted photon density, and $dk_i$ is the differential solid angle.

$\rho$ is constant for diffuse surfaces.

### Surface reflectance types

- diffuse (D): light scatters in all directions equally
- rough specular (glossy) (S): "specular lobe", _mostly_ bouncing in the mirror direction
- mirror (ideal specular) (S): Exact mirror reflections

Heckbert light path notation:
    - L (light)
    - D (diffuse surface)
    - S (specular/glossy/mirror surface)
    - E (eye or camera)

Regular expressions:
    - A|B (A or B)
    - A* (zero or more A's)
    - A+ (one or more A's)

One bounce model: Path = LDE

Multiple bounces (closed room, diffuse environment): LD*E

Classical: LDE is calculated, LSDE is omitted (how to fix?)

Caustics: Bright spot from bending through lens.
- Path: LSSDE (SS = into and out of lens)
- Classical ray tracing has no caustics


## Mon Mar 11

### The rendering equation

Jim Kajiya, Giggraph 1986

Given:
- scene geometry
- BRDF's
- light sources (emitting geometry)

Determine:
- light that reaches the eye

Approach: start from the refl. equation (integral over directions), and convert to integral over points in scene
- solid angles -> areas on surface
- add visibility term $V(x, x') \in (0, 1)$ (is the path blocked?)
- add light sources
- some directions turn into points

$$dk_i = \frac{\cos{\theta_o} dA}{|x-x'|^2}$$
$G(x,x') = \frac{\cos{\theta_i}\cos{\theta_o}}{|x-x'|^2}$$
"geometry" term


$x'$ is the surface reflection point, $x$ is the point in the scene, $dA$ is the differential surface area, $\theta_o$ is the surface normal angle

$$L(x',k') = L_e(x',k') + \int_{x \in S}{\rho(k,k')L(x,k)V(x,x')G(x,x')dA}$$

- $L$: outgoing light from $x'$ to eye
- $L_e$: emitted light from $x'$
- $x \in S$: over all surfaces
- $\rho$: reflectence (brdf)
- $V$: visibility term
- $G$: geometry term
- $dA$: differential area

```
Radiance get_radiance(ray):
    cast ray into scene, hit point x with direction k
    radiance = emit_radiance(x, k)
    for each light source x' in direction k':
        if sees_light(x, light):
            radiance += weight_radiance(brdf(k,k'), light_radiance)

    if random < reflect_prob(x):
        create reflected_ray, direction k'
        radiance += weight_radiance(brdf(k, k'), get_radiance(reflected_ray))
    else
        create refracted_ray, direction k'
        radiance += weight_radiance(brdf(k, k'), get_radiance(refracted_ray))

    return radiance
```

## Wed Mar 13

### Photon mapping

#### Pass 1

- trace photons from light
- store photons at surfaces
- effects:
    - caustics
    - indirect diffuse illumination
    - participating media (e.g. shafts of light in fog)

Caustics: light focussed by lens/mirror (specular surfaces)
1) shoot photons from lights
2) view scene from eye, show concentrations of photons

Photons:
- amount of energy
- store location where they hit on diffisue surface
- store direction of travel (when we have non-diffuse surfaces)

Store photons using kd-trees. works for _any_ shaped surfaces

##### k-D trees

Support fast proximity queries (e.g. gimme the 10 nearest photons)

- stored split planes
- cycle through x,y,z
- stop splitting when last split has few photons

#### Caustics Pass 2

- shoot rays from eye
- calcuate diffuse illum. normally ($c = c_r c_l (N \cdot L)$)
- add the photons' contribution - photon density:
    - find nearest N photons (e.g. N=50)
    - density proportional to needed search radius
- use thise to determine if it's a bright spot

Indirect diffuse Illumination
- Pass 1
1) shoot photons from light, bounce them off diffuse survaces (in a random direction)
    - don't store the first hit from the light source
    - photons can bounce more than once
    - bounce prob. based on reflectance of material (Russian Roulette sampling)
2) render from eye, adding indirect term
- Pass 2
- Final gather
    - k_eye - density estimation at each scene point to estimate L(x_1,k_1) (light leaving)
    - Irradiance: $L(x,k_{eye}) = \sum_i{\rho(k_i,k_{eye})L(x,k_i)}$
    - Irradiance caching (Greg Ward) - need more detail at e.g. wall corners, under objects, so only do final gather in those places and interpolate
        - Estimate distance + normal change at nearby cache points


## Fri Mar 15

Anti-Aliasing

Aliasing: High frequency components of a sampled signal misidentified as lower frequencies.

Jagged edges of polygons is an example of _spatial_ aliasing.
Also moire patterns of grids.
- samples: pixel colors
- signal: ideal image
- frequency components: regularly spaced pixels
    - isn't the analagous to the _time domain_, not the _frequency domain_?
      Frequency domain analog would be the color it's aliased to.
Corrected:
- time components: pixel spacing
- frequency component: pixel color

- use grayscale to correct
- box filter (not ideal)
    - area-weighted filter
    - value = volume under the box

triangle filter give us sub-pixel precision


## Mon Mar 25

Project, rasterize, z=buffer

Perspective projection:

y'/1 = y/|z|

x' = x/|z|

y' goes up to p', the height y' where we see point p intersecting the view plane at y' = p'


rasterization: turn simple geometry into pixels
(e.g. lines, triangles, ...)
- given: 2D  vertices of triangle
- calc pixels in triangle, and find min+max x & y coords
- then, fill in each pixel in the area _inside the triangle_


z-buffer: visible surface determination

framebuffer: pixel colors (2d array)
z-buffer: z or depth value per-pixel

as we rasterize, each pixel has "depth" recorded

for each pixel in window (x, y)
    setz(x, y, far)

for each triangle in scene
    for each pixel in triangle
        if pz < bufferz(x,y)
            write_pixel(x,y, r,g,b)
            setz(x,y, pz)

ray tracing:
for each pixel x,y
    create ray R
    for each object o_i
        is o_i closest along ray R?

rasterization
for each triangle o_i
    for each pixel in triangle o_i
        is pz closest at x,y?


## Wed Mar 27

Vertex processors (programmable!)
- can calculate per-vertex shading

Rasterizer: (not programmable)
- converts screen-space geometry (e.g. triangles) to fragments
- interpolates per-vertex info


NVidia GPUs

- gpu chip
    - (many streaming multiprocessors, executing different instructions)
        - each multiprocessors has many cores ("CUDA cores")
        - SIMD (each core acts on different vertices/fragments but doing the same operations)


## Mon Apr 1

Rexture maps

texture - vary color across surface

examples:
- wood grain
- photo on book jacket
- carpet pattern

texture space / screen space

texture (image) space: s/t (or u/v), s \in [0,1], t \in [0,1] (for example)

texel (short for "texture element") ~= pixel

Render textured polygons
1) Interpolate s and t across polygon
2) Find color from texture at (s,t) - texture lookup
3) Color the given pixel using texel color and shading

Interpolate s and t across poly: barycentric coordinates

triangle p1, p2, p3, with barycentric coords (s1,t1), (s2,t2), (s3,t3)
p in center. \alpha, \beta, \gamma are the areas of the three sections of the triangle split by p

s = \alpha s1 + \beta s2 + \gamma s3

Perspective: x' = x/|z|, y' = y/|z|, z' = 1/|z|

1) divide s and t by z to give s' and t'
2) interpolate s' and t' across poly. (barycentric)
3) divide s' and t' by z' (per pixel) -> s and t
4) perform texture lookup and shade

if we round texture coords to integer: locations = nearest neighbor lookup
- gives blocky/pixelated results

instead, average over the 4 nearest texels to given point
- bilinear/trilinear interpolation: interp across x dim, then y
    r = r0 + b(r1 - r0), with r1 = r01 + a(r11 - r01), and r0 = r00 + a(r10 - r00)

Texture minification - avoid sparkling, scintillation
- create multiple textures at different resolutions, low pass filter
- mipmap - image pyramid + trilinear interpolation (bilinear across texture dims and then linear across neighbor in image hierarchy)
- how to pick resolutions?
    - want: one pixel on screen space
    - square pixel in screen space maps to a quadrilateral in texture space.
    - pixel footprint: turn the quadrilateral projection into an elipse, calculate a differential step length d, and choose texture size based on d


## Wed Apr 3

volume data: values given at each 3d cell, "sampled" 3d data
typical grid sizes: 512x512, 128x128

Volume data sources:
- Medical CT scan - permiability to x rays
- MRI - measure amount of hydrogen in tissue
- Fluid simulation - pressure, vorticity, speed
- Engineering simulation - temp, stress, strain
- Smoke, clouds

Volume rendering
- isosurface extraction (e.g. marching cubes)
- direct volume rendering (ray tracing)

Direct volume rendering
1) classify each voxel (gives reflectance, opacity (\alpha, prob. light will hit a particle inside))
2) calculate shading (normals + reflectance -> radiance r)
3) render the voxels

voxel composition equation

final color = $\alpha_1 r_1 + (1-\alpha_1)[\alpha_2 r_2 + (1-\alpha_2)[\alpha_3 r_3 + (1-\alpha_3)\cdots]]$


classification: voxel scalar value -> opacity, reflectance

surface normals are determined from gradient values
- gradient is perp. to iso-contours
- $G$ = gradient, $N = -G/\|G\|$
- use finite differences to estimate gradients
- $G = \partial v/dx, \partial v/dy$


Edge Collapse
- repeat until N vertices remain
- select edge with lowest cost (priority queue)
- collapse edge
- pick new vertex position
ce new vertex?
- midpoint of edge
- one of the o.g. verts
- place causing the least surface movement

removes: 1 vert, 2 tris, 3 edges

cost of edge collapse: dist. to plances of nearby tri's

signed distance to line
- f(x,y) = ax + by + c = 0

sq. dist to line
- g(x,y) = (ax + by + c)^2 = 0

g+H

line l = ax + by + c = L^TX

(L^TX)^2 = (ax + by + c)^2

Plane P = [a b c d], x = [x y z 1]

signed function f(x) = ax + by + cz + d = P^TX

squared dist: g(x) = (P^TX)^2 = X^T(PP^T)x (4x4 matrix, K_P = PP^T is symmetric, distance sq. to _any_ number of planes)

error: g(x) = \sum_{i=1}^n{X^T K_P_i X} = X^T (\sum K_p_i) X (\sum K_p_i \trangleeq M is 4x4 matrix)

g(x) = X^TMX

minimize x for given M: g(x) = 1/2 X^TAX + X^T + c, (A is 3x3 mat., x^T is a 3-vector, c is scalar)

\Delta g = Ax - b = 0 -> X = A^{-1}b



## Wed Apr 10

Voronoi diagrams

given: set of points (sites) p1,p2,...pn in the plane
create: partition of the plane into regions by assigning each point to its nearest site

Fast approximation: rasterize cones


even point placement in 2d: repeat
1) place points at random in plane
2) calc Voronoi regions (using cones & gpu)
3) find centroids of regions
4) move each point to its region's centroid

calculate Voronoi diagrams exactly: several methods, fastes iones are O(nlogn)

Fortune's Alg:
- sweep a horizontal line from top to bottom
- maingain 1d slice thru 2d surface

## Fri Apr 12

Delauney triangularion

Straight line dual of the Voronoi Diagram

Connect the _adjacent_ sites from Voronoi diagram

Delauney property: Circumcircles of triangle vertices contain no other sites
- Can use this prop. to calculate triangles via edge-flips
    - When a side is found inside a circumcircle of a triangle, flip the edge shared with an adjascent triangle
    - After flippin middle edge, neighter triangle will have points in the other's circumcircle
    - Triangles will be less skinny (closer to equilateral)

Calc Delauney triangularion by adding one point at a time
1) enclose points in 2 triangles
2) add one point to triangulation
    - find which triangle contains the point
    - place inside triangle, make 3 new triangles (throw "parent" out)
    - see if nearby edges need to be flipped

Repeaat (2) until all points added

DT has nice property:
- has the largest minimum angle of all possible triangulations of the sides

DT uses:
- finite element analysis
- terrain visualization
- graph algorithms

More useful with _constraint edges_: Constrained DT (CDT)
- Given: Set of sites p1, p2, ..., pn, and adges e1, e2, ..., em
- Make: Triangulation of sites that includes all constraint edges + maintains Delaunay property when possible.

Algo:
1) calc. DT w/ no constraints
2) repeatedly add constraint edges
    - delete triangles on the edge's path
    - insert constraint edge
    - re-triangulate the two resulting holes
