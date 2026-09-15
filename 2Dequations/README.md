MATLAB code for approximating solutions of Lu=f(u), where L is an elliptic differential operator. 

The methods approximate critical points of the functional
S(u)=(1/2)*int Lu*u dx - int F(u) dx, where F'=f.

All three methods seek local minimizers of S on the Nehari manifold 
N={u : P(u)=0}
where P(u)=<S'(u),u>.

gradsolver is a constrained gradient descent method

rgradsolver is a constrained Rimennan gradient descent method

rnewtsolver is a constrained Riemannian inexact Newton method, where the inverse of the Hessian of S is approximated by a truncated Neumann series.


