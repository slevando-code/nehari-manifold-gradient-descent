function E=gradsolver_diff(p,q,b,c,A,Nx,tol,u,delta,plot_final)

% Nehari manifold gradient descent method for approximating solutions of
% u''''+bu''+cu=|u|^(q-1)*u-|u|^(p-1)*u where q>p
% Sample usage E=gradsolver_diff(3,5,1.8,1,50,1000,1e-10,0,.5,1)
% Spatial interval is [-A,A], with Nx subintervals.
% delta = step size in gradient descent step
% u = initial guess, set u=0 to use default Gaussian
% tol = error tolerance for H^2 norm of gradient of S
% Output E is vector of H^2 errors between successive iterates
% set plot_final==1 to plot final iterate
% c must be positive and b must be less than 2*sqrt(c)

arguments
    p (1,1) double = 3
    q (1,1) double = 5
    b (1,1) double = 1.8
    c (1,1) double = 1
    A (1,1) double = 50
    Nx (1,1) double {mustBeInteger, mustBePositive} = 1000
    tol (1,1) double {mustBePositive} = 1e-10
    u (1,:) double = 0
    delta (1,1) double {mustBePositive} = 0.5
    plot_final double = 1
end

if c<=0 || b>=2*sqrt(c)
    error('Choose c>0 and b<2*sqrt(c).')
end

if p>=q
    error('Choose p<q.')
end

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/Nx)*(2*(1:Nx)-Nx); % Real spatial variable

E=[];

if u==0
    u=exp(-x.^2);
end

m = pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;

% Determine alpha for which P(alpha*u)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
N1=sum((abs(u).^(q+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=PRoot_diff(p,q,2*I,N2,N1);

u=alpha*u;

err=tol+1;

while (err>tol)

    Oldu=u;

    % Compute gradient of S at u
    f=abs(u).^(q-1).*u-abs(u).^(p-1).*u;
    gradS=real(u-ifft(fft(f)./m));

    % Step in direction -gradS
    u=u-delta*gradS;

    % Determine alpha for which P(alpha*u)=0
    N2=sum((abs(u).^(p+1)))*(A/Nx);
    N1=sum((abs(u).^(q+1)))*(A/Nx);
    dI=real(ifft(m.*fft(u)));
    I=(1/2)*sum(dI.*u)*(A/Nx);
    alpha=PRoot_diff(p,q,2*I,N2,N1);
    u=alpha*u;

    % Compute H^2 norm of gradient of S
    LgradS=real(ifft(m.*fft(gradS)));
    err=sqrt(sum(LgradS.*gradS)*(A/Nx));

    % Compute H^2 norm of difference between iterates
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


%Newton's Method for approximating the positive root of f(x)=a+b*x^(p-1)-c*x^(q-1)
function x=PRoot_diff(p,q,a,b,c)
tol=1e-14;
x=(b/c)^(1/(q-p));
e=1;

while (abs(e)>tol)
    f=a+b*x^(p-1)-c*x^(q-1);
    df=b*(p-1)*x^(p-2)-c*(q-1)*x^(q-2);
    e=f/df;
    x=x-e;
end
