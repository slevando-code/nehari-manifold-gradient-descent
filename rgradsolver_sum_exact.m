function E=rgradsolver_sum_exact(p,b,A,Nx,tol,u,delta,plot_final,plot_error)

% Riemannian gradient descent method for approximating solutions of
% u''''+bu''+cu=|u|^(p-1)*u+|u|^(q-1)*u
% Sample usage: E=rgradsolver_sum_exact(3,-4,100,4000,1e-10,0,1.1,1,1)
% Spatial interval is [-A,A], with Nx subintervals.
% delta = step size in gradient descent step
% u = initial guess, set u=0 to use default Gaussian
% tol = H^2 norm error tolerance
% output E is vector of H^2 errors between iterates and exact solution
% set plot_final==1 to plot final iterate
% set plot_error==1 to plot final difference
% b must be less than -sqrt(p*(3*p-1)/(p+1))

if b >= -sqrt(p*(3*p-1)/(p+1))
    error('b must be less than -sqrt(p*(3*p-1)/(p+1))')
end

E=[];

q=2*p-1;

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/(Nx))*(2*(1:Nx)-Nx); % Real spatial variable

% Compute exact solution
r=2/(p-1);
d=sqrt((-sqrt(r*(r+1)*(r+2)*(r+3))-b*r*(r+1))/(2*r*(r+1)*(r^2+2*r+2)));
M=(sqrt(r*(r+1)*(r+2)*(r+3))*d^2)^(1/(p-1));
c=-r^4*d^4-b*r^2*d^2;
exact=M*sech(d*x).^r;


if u==0
    u=exp(-x.^2);
end

v= pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;

% Determine alpha for which P(alpha*Phi)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
N1=sum((abs(u).^(q+1)))*(A/Nx);
dI=real(ifft(v.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=((-N2+sqrt(N2^2+8*I*N1))/(2*N1))^(1/(p-1));

u=alpha*u;

H2err=tol+1;

while (H2err>tol)


    Oldu=u;

    %Compute gradient of S at Phi
    f=abs(u).^(q-1).*u+abs(u).^(p-1).*u;
    gradS=real(u-ifft(fft(f)./v));

    %Compute gradient of P at Phi
    df = q*abs(u).^(q-1)+p*abs(u).^(p-1);
    gradP=real(2*u-ifft(fft(f+u.*df)./v));

    %Compute projection of gradient of S onto tangent space of Nehari manifold
    %(Riemannian Gradient)
    proj = gradS-real(sum(v.*fft(gradS).*conj(fft(gradP)))/sum(v.*fft(gradP).*conj(fft(gradP))))*gradP;

    %Experiment: Try projecting onto Phi perp
    %proj = gradS-real(sum(v.*fft(gradS).*conj(fft(Phi)))/sum(v.*fft(Phi).*conj(fft(Phi))))*Phi;


    % Move in direction of Riemannian gradient
    u=u-delta*proj;

    N2=sum((abs(u).^(p+1)))*(A/Nx);
    N1=sum((abs(u).^(q+1)))*(A/Nx);
    dI=real(ifft(v.*fft(u)));
    I=(1/2)*sum(dI.*u)*(A/Nx);

    % Determine alpha for which P(alpha*Phi)=0
    alpha=((-N2+sqrt(N2^2+8*I*N1))/(2*N1))^(1/(p-1));
    u=alpha*u;

    diff=u-exact;
    Ldiff=real(ifft(v.*fft(diff)));
    H2err=sqrt(sum(Ldiff.*diff)*(A/Nx))
    E=[E ; H2err];

end

if plot_final==1
    figure(1)
    hold off
    plot(x,u,'b')
    hold on
    plot(x,exact,'r')
end

if plot_error==1
    figure(2)
    hold off
    plot(x,u-exact)
end


