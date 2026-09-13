function E=rnewtsolver_sum_exact(p,b,A,Nx,tol,u,delta,n,plot_final,plot_error)

% Nehari manifold inexact Newton method for approximating solutions of
% u''''+bu''+cu=|u|^(p-1)*u+|u|^(q-1)*u where q=2*p-1
% Sample usage: E=rnewtsolver_sum_exact(3,-4,100,4000,1e-10,0,1.2,2,1,1)
% Spatial interval is [-A,A], with Nx subintervals.
% delta = step size in approximation of Hessian of S
% u = initial guess, set u=0 to use default Gaussian
% tol = H^2 norm error tolerance
% output E is vector of H^2 errors between iterates and exact solution
% n = # of terms in Neumann series approximation of inverse of Hessian of S
% set plot_final==1 to plot final iterate
% set plot_error==1 to plot final difference
% b must be less than -sqrt(p*(3*p-1)/(p+1))

if b >= -sqrt(p*(3*p-1)/(p+1))
    error('b must be less than -sqrt(p*(3*p-1)/(p+1))')
end

q=2*p-1;

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/(Nx))*(2*(1:Nx)-Nx); % Real spatial variable


% Calculate c and parameters in exact solution in terms of b and p
r=2/(p-1);
d=sqrt((-sqrt(r*(r+1)*(r+2)*(r+3))-b*r*(r+1))/(2*r*(r+1)*(r^2+2*r+2)));
M=(sqrt(r*(r+1)*(r+2)*(r+3))*d^2)^(1/(p-1));
c=-r^4*d^4-b*r^2*d^2;


% Exact solution
exact=M*sech(d*x).^r;

if u==0
    u=exp(-x.^2);
end

m = pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;


% Determine alpha for which P(alpha*u)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
N1=sum((abs(u).^(q+1)))*(A/Nx);
dJ=real(ifft(m.*fft(u)));
J=(1/2)*sum(dJ.*u)*(A/Nx);
alpha=((-N2+sqrt(N2^2+8*J*N1))/(2*N1))^(1/(p-1));

u=alpha*u;


% Initialize counter and errors

H2err=tol+1;
E=[];

while (H2err>tol)


    %Compute gradient of S at u
    f=abs(u).^(q-1).*u+abs(u).^(p-1).*u;
    gradS=real(u-ifft(fft(f)./m));

    %Compute gradient of P at u
    df = q*abs(u).^(q-1)+p*abs(u).^(p-1);
    gradP=real(2*u-ifft(fft(f+u.*df)./m));



    %Project onto tangent space of N

    w0 = gradS-real(sum(m.*fft(gradS).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
    w=w0;

    for i=1:n
        J=ifft(fft(df.*w)./m);
        J=J-real(sum(m.*fft(J).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
        w=w0+(1-delta)*w+delta*J;
    end


    u=u-delta*w;



    N2=sum((abs(u).^(p+1)))*(A/Nx);
    N1=sum((abs(u).^(q+1)))*(A/Nx);
    dI=real(ifft(m.*fft(u)));
    I=(1/2)*sum(dI.*u)*(A/Nx);

    % Determine alpha for which P(alpha*u)=0
    %alpha=PRoot(p,q,2*I,N2,N1);
    alpha=((-N2+sqrt(N2^2+8*I*N1))/(2*N1))^(1/(p-1));

    u=alpha*u;

    diff=u-exact;
    Ldiff=real(ifft(m.*fft(diff)));
    H2err=sqrt(sum(Ldiff.*diff)*(A/Nx));
    E=[E ; H2err];

end


% Plot exact solution together with final iterate
if plot_final==1
    figure(1)
    hold off
    plot(x,u,'b','LineWidth',1)
    hold on
    plot(x,exact,'r')
end




%Plot exact error of last iterate
if plot_error==1
    figure(2)
    hold off
    plot(x,u-exact,'k')
end

