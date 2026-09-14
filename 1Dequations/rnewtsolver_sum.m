function E=rnewtsolver_sum(p,q,b,c,A,Nx,tol,u,delta,n,plot_final)

% Nehari manifold inexact Newton method for approximating solutions of
% u''''+bu''+cu=|u|^(p-1)*u+|u|^(q-1)*u
% Sample usage: E=rnewtsolver_sum(3,6,2,2,50,1000,1e-10,0,1.2,2,1)
% Spatial interval is [-A,A], with Nx endpoints. Nx must be even.
% delta = step size in approximation of Hessian of S
% u = initial guess. Set u=0 to use default Gaussian
% tol = error tolerance for H^2 norm of gradient of S
% n = # of terms in Neumann series approximation of inverse of Hessian of S
% Output E is vector of H^2 errors between successive iterates
% set plot_final==1 to plot final iterate
% c must be positive and b must be less than 2*sqrt(c)

if c<=0 || b>=2*sqrt(c)
    error('Choose c>0 and b<2*sqrt(c)')
end

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/(Nx))*(2*(1:Nx)-Nx); % Real spatial variable

% If no initial guess specified, use this
if u==0
    u=exp(-x.^2);
end

% Fourier multiplier that corresponds to linear operator Lu=u''''+bu''+cu
m= pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;

% Determine alpha for which P(alpha*Phi)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
N1=sum((abs(u).^(q+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=PRoot_sum(p,q,2*I,N2,N1);
u=alpha*u;

err=tol+1;
E=[];

while (err>tol)

    Oldu=u;

    f=abs(u).^(q-1).*u+abs(u).^(p-1).*u;
    df=p*abs(u).^(p-1)+q*abs(u).^(q-1);

    gradS=real(u-ifft(fft(f)./m));

    %Compute gradient of P at u
    df=p*abs(u).^(p-1)+q*abs(u).^(q-1);
    gradP=real(2*u-ifft(fft(f+u.*df)./m));

    %Project onto tangent space of Nehari manifold

    w0 = gradS-real(sum(m.*fft(gradS).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
    w=w0;

    % Approximate inverse of Hessian of S at u by truncated Neumann series

    for i=1:n
        J=ifft(fft(df.*w)./m);
        J=J-real(sum(m.*fft(J).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
        w=w0+(1-delta)*w+delta*J;
    end

    u=u-delta*w;

    % Determine alpha for which P(alpha*Phi)=0
    N2=sum((abs(u).^(p+1)))*(A/Nx);
    N1=sum((abs(u).^(q+1)))*(A/Nx);
    dI=real(ifft(m.*fft(u)));
    I=(1/2)*sum(dI.*u)*(A/Nx);
    alpha=PRoot_sum(p,q,2*I,N2,N1);
    u=alpha*u;

    LgradS=real(ifft(m.*fft(gradS)));
    err=sqrt(sum(LgradS.*gradS)*(A/Nx));

    diff=u-Oldu;
    Ldiff=real(ifft(m.*fft(diff)));
    differr=sqrt(sum(Ldiff.*diff)*(A/Nx));
    E=[E;differr];

end

% Plot final iterate
if plot_final==1
    figure(1)
    hold off
    plot(x,u,'LineWidth',1)
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
    ylabel({'$u$'},'interpreter','latex','FontSize',12)
end

function x=PRoot_sum(p,q,a,b,c)
%Newton's Method for approximating the positive root of f(x)=a-b*x^(p-1)-c*x^(q-1)
tol=1e-15;
x=(a/b)^(1/(p-1));
e=1;

while (abs(e)>tol)
    f=a-b*x^(p-1)-c*x^(q-1);
    df=-b*(p-1)*x^(p-2)-c*(q-1)*x^(q-2);
    e=f/df;
    x=x-e;
end
