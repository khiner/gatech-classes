# Graphics notes

## Fri Jan 12

### Eye rays

x(t) = x_0 + t(x_1 - x_0)

Translate x \in [0,w] -> x' \in [-k,k] (the view plane for hw 1A)

x' = (x - w/2) \* (2k/w)
y'

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
