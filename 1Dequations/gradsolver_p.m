function E=gradsolver_p(p,b,c,A,Nx,tol,u,delta,plot_final)

% Nehari gradient descent method for approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Sample usage:
% E=gradsolver_p(3,1.7,1,50,1000,1e-10,0,1.2,1)
% Spatial interval is [-A,A], with Nx subintervals. Nx must be even.
% delta = step size in gradient descent step
% u = initial guess. Set u=0 to use default Gaussian
% tol = error tolerance for H^2 norm of gradient of S
% Output E is vector of H^2 norms of differences between iterates
% c must be positive and b must be less than 2*sqrt(c)
% set plot_final==1 to plot final iterate

if c<=0 || b>=2*sqrt(c)
    error('Choose c>0 and b<2*sqrt(c)')
end

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/Nx)*(2*(1:Nx)-Nx); % Real spatial variable

if u==0
    u=exp(-x.^2);
end

% Fourier multiplier that corresponds to linear operator Lu=u''''+bu''+cu
m= pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;

% Determine alpha for which P(alpha*Phi)=0
N=sum((abs(u).^(p+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=(2*I/N)^(1/(p-1));
u=alpha*u;

E=[];
Reserr=tol+1;

while (Reserr>tol)

    Oldu=u;

    f=abs(u).^(p-1).*u;

    gradS=real(u-ifft(fft(f)./m));

    % Step in direction of gradient
    u=u-delta*gradS;

    % Determine alpha for which P(alpha*Phi)=0 and rescale
    N=sum((abs(u).^(p+1)))*(A/Nx);
    dI=real(ifft(m.*fft(u)));
    I=(1/2)*sum(dI.*u)*(A/Nx);
    alpha=(2*I/N)^(1/(p-1));
    u=alpha*u;

    % Compute H^2 norms of gradient of S and difference between u and Oldu
    LgradS=real(ifft(m.*fft(gradS)));
    Reserr=sqrt(sum(LgradS.*gradS)*(A/Nx));

    diff=u-Oldu;
    Ldiff=real(ifft(m.*fft(diff)));
    Differr=sqrt(sum(Ldiff.*diff)*(A/Nx));

    E=[E;Differr];

end

% Plot final iterate
if plot_final==1
    figure(1)
    hold off
    plot(x,u,'LineWidth',1)
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
    ylabel({'$u$'},'interpreter','latex','FontSize',12)
end


