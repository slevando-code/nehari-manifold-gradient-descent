function E=rnewtsolver_p_exact(p,b,A,Nx,tol,u,delta,n,plot_final,plot_error)

% Inexact Newton method for approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Sample usage: 
% E=rnewtsolver_p_exact(3,-2.5,50,2000,1e-10,0,1.2,2,1,1)
% Spatial interval is [-A,A], with Nx subintervals. Nx must be even.
% delta = step size in gradient descent step
% u = initial guess
% tol = H^2 norm error tolerance between u and exact solution
% n = # of terms in Neumann series approximation of inverse of Hessian of S
% b must be less than 0
% Output E is vector of H^2 errors between iterates and exact solution
% set plot_final==1 to plot final iterate
% set plot_error==1 to plot final difference

if b >= 0
    error('b must be less than zero')
end


% Calculate c and parameters in exact solution in terms of b and p
c=(2*(p+1)*b/(p^2+2*p+5))^2;
d=((p-1)*c^(1/4))/sqrt(8*(p+1));
M=(((p+3)*(3*p+1)*c)/(8*(p+1)))^(1/(p-1));
r=4/(p-1);

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/Nx)*(2*(1:Nx)-Nx); % Real spatial variable

% Exact solution
exact=M*sech(d*x).^r;


if u==0
 u=exp(-x.^2);
end


m = pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;



% Determine alpha for which P(alpha*Phi)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=(2*I/N2)^(1/(p-1));
u=alpha*u;

% Initialize counter and errors

H2err=1;
E=[];

while (H2err>tol)

f=abs(u).^(p-1).*u;
df=p*abs(u).^(p-1);


% Compute gradient of S at u
gradS=real(u-ifft(fft(f)./m));

% Compute gradient of P at u
gradP=real(2*u-(p+1)*ifft(fft(f)./m));

% Project gradS onto tangent space of N at u 
w0 = gradS-real(sum(m.*fft(gradS).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
w=w0;

% Approximate inverse of Hessian of S at u by truncated Neumann series

for i=1:n
   J=ifft(fft(df.*w)./m);
   J=J-real(sum(m.*fft(J).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
   w=w0+(1-delta)*w+delta*J;
end

% Take approximate Newton iteration step.
u=u-delta*w;

% Determine alpha for which P(alpha*u)=0 and rescale u
N2=sum((abs(u).^(p+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=(2*I/N2)^(1/(p-1));
u=alpha*u;

Exactdiff=u-exact;
LExactdiff=real(ifft(fft(Exactdiff).*m));
H2err=sqrt(sum(LExactdiff.*Exactdiff)*(A/Nx));

E=[E; H2err];
end

% Plot exact solution together with final iterate
if plot_final==1
    figure(1)
    hold off
    plot(x,u,'b','LineWidth',1)
    hold on
    plot(x,exact,'r')
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
    ylabel({'$u$'},'interpreter','latex','FontSize',12)
end


%Plot exact error of last iterate
if plot_error==1
    figure(2)
    plot(x,u-exact,'k')
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
end

