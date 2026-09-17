# nehari-manifold-gradient-descent
MATLAB code for solving elliptic equations using Nehari manifold gradient descent algorithm.

A collection of m-files for approximating nontrivial solutions of equations of the form Lu=f(u), where L is an elliptic operator.

Solutions are critical points of the functional S defined by
S(u)=(1/2) int_R^n Lu * u dx - int_R^n F(u) dx
where F'=f. 

Thus all solutions lie on the Nehari manifold P(u)=0, where P(u)=<S'(u),u>. There are three methods, all of which involve minimizing S constrained to the Nehari manifold. 

The first method, "gradsolver", is a simple gradient descent followed by a scaling back to the Nehari manifold.
The second method, "rgradsolver", is a Riemannian gradient descent, where the gradient descent direction is determined by projecting the gradient of S onto the tangent space of the Nehari manifold. 
The third method, "rnewtsolver", is an inexact Newton method in which the Newton iteration step is approximated by using a truncated Neumann series to approximate the Hessian of S.

For each method:

_p stands for the pure power nonlinearity f(u)=|u|^(p-1)*u

_sum stands for a sum of powers f(u)=|u|^(p-1)*u+|u|^(q-1)*u

_diff stands for a difference of powers f(u)=|u|^(q-1)*u-|u|^(p-1)*u where p<q





