function E=gradsolver_sum_exact(p,b,A,Nx,tol,u,delta,plot_final,plot_error)

% Gradient descent method for approximating solutions of
% u''''+bu''+cu=|u|^(p-1)*u+|u|^(q-1)*u where q=2*p-1
% Sample usage: E=gradsolver_sum_exact(3,-4,100,1000,1e-10,0,1.2,1,1)
% Spatial interval is [-A,A], with Nx subintervals
% delta = step size in gradient descent step
% u = initial guess, set u=0 to use default Gaussian
% tol = H^2 norm error tolerance between u and exact solution
% set plot_final==1 to plot final iterate
% set plot_error==1 to plot final difference between u and exact solution
% b must be less than -sqrt(p*(3*p-1)/(p+1))
% Output E is vector of H^2 errors between successive iterates

if b >= -sqrt(p*(3*p-1)/(p+1))
    error('b must be less than -sqrt(p*(3*p-1)/(p+1))')
end

q=2*p-1;

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/(Nx))*(2*(1:Nx)-Nx); % Real spatial variable

% Exact solution
r=2/(p-1);
d=sqrt((-sqrt(r*(r+1)*(r+2)*(r+3))-b*r*(r+1))/(2*r*(r+1)*(r^2+2*r+2)));
M=(sqrt(r*(r+1)*(r+2)*(r+3))*d^2)^(1/(p-1));
c=-r^4*d^4-b*r^2*d^2;

exact=M*sech(d*x).^r;

if u==0
    u=exp(-x.^2);
end

m= pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;

% Determine alpha for which P(alpha*Phi)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
N1=sum((abs(u).^(q+1)))*(A/Nx);
dJ=real(ifft(m.*fft(u)));
J=(1/2)*sum(dJ.*u)*(A/Nx);
alpha=((-N2+sqrt(N2^2+8*J*N1))/(2*N1))^(1/(p-1));

u=alpha*u;

H2err=tol+1;

E=[];


while (H2err>tol)

    %Compute gradient of S at Phi
    f=abs(u).^(q-1).*u+abs(u).^(p-1).*u;
    gradS=real(u-ifft(fft(f)./m));

    % Move in direction of Riemannian gradient
    u=u-delta*gradS;

    % Determine alpha for which P(alpha*Phi)=0
    N2=sum((abs(u).^(p+1)))*(A/Nx);
    N1=sum((abs(u).^(q+1)))*(A/Nx);
    dJ=real(ifft(m.*fft(u)));
    J=(1/2)*sum(dJ.*u)*(A/Nx);
    alpha=((-N2+sqrt(N2^2+8*J*N1))/(2*N1))^(1/(p-1));

    u=alpha*u;

    diff=u-exact;
    Ldiff=real(ifft(m.*fft(diff)));
    H2err=sqrt(sum(Ldiff.*diff)*(A/Nx));
    E=[E;H2err];

end

if plot_final==1
    figure(1)
    hold off
    plot(x,u,'b','LineWidth',1)
    hold on
    plot(x,exact,'r')
end

if plot_error==1
    figure(2)
    hold off
    plot(x,u-exact,'LineWidth',1)
end

