function E=rgradsolver_sum(p,q,b,c,A,Nx,tol,u,delta,plot_final)

% Nehari manifold Riemannian gradient descent method for approximating solutions of
% u''''+bu''+cu=|u|^(p-1)*u+|u|^(q-1)*u
% Sample usage:
% E=rgradsolver_sum(3,6,2,2,50,1000,1e-10,0,1.2,1)
% Spatial interval is [-A,A], with Nx endpoints. Nx must be even.
% delta = step size in gradient descent step
% u = initial guess, enter 0 to use Gaussian
% tol = error tolerance for H^2 norm of gradient of S
% Output E is vector of H^2 errors between succesive iterates
% set plot_final==1 to plot final iterate
% c must be positive and b must be less than 2*sqrt(c)

if c<=0 || b>=2*sqrt(c)
    error('Choose c>0 and b<2*sqrt(c)')
end

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/(Nx))*(2*(1:Nx)-Nx); % Real spatial variable

if u==0
    u=exp(-x.^2);
end

m = pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;

% Determine alpha for which P(alpha*u)=0 and scale u
N2=sum((abs(u).^(p+1)))*(A/Nx);
N1=sum((abs(u).^(q+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=PRoot_sum(p,q,2*I,N2,N1);
u=alpha*u;

E=[];
H2res=tol+1;

while (H2res>tol)

    Oldu=u;

    %Compute gradient of S at u
    f=abs(u).^(q-1).*u+abs(u).^(p-1).*u;
    gradS=real(u-ifft(fft(f)./m));

    %Compute gradient of P at u
    df = q*abs(u).^(q-1)+p*abs(u).^(p-1);
    gradP=real(2*u-ifft(fft(f+u.*df)./m));

    %Compute projection of gradient of S onto tangent space of Nehari manifold
    %(Riemannian Gradient)
    proj = gradS-real(sum(m.*fft(gradS).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;

    % Move in direction of Riemannian gradient
    u=u-delta*proj;


    % Determine alpha for which P(alpha*Phi)=0
    N2=sum((abs(u).^(p+1)))*(A/Nx);
    N1=sum((abs(u).^(q+1)))*(A/Nx);
    dI=real(ifft(m.*fft(u)));
    I=(1/2)*sum(dI.*u)*(A/Nx);
    alpha=PRoot_sum(p,q,2*I,N2,N1);
    u=alpha*u;

    LgradS=real(ifft(m.*fft(gradS)));
    H2res=sqrt(sum(LgradS.*gradS)*(A/Nx));
    
    diff=Oldu-u;
    Ldiff=real(ifft(m.*fft(diff)));
    H2diff=sqrt(sum(Ldiff.*diff)*(A/Nx));
    E=[E; H2diff];
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
