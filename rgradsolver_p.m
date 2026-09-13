function EDiff=rgradsolver_p(p,b,c,A,Nx,tol,u,delta,plot_final)

% Riemannian gradient descent method for approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Sample usage E=rgradsolver_p(3,1.8,1,50,1000,1e-10,0,1.2,1)
% Spatial interval is [-A,A], with Nx subintervals. Nx must be even
% delta = step size in gradient descent step
% u = initial guess, set u=0 to use default Gaussian
% tol = error tolerance for H^2 norm of gradient of S
% Output E is vector of H^2 norms of differences between iterates
% c must be positive and b must be less than 2*sqrt(c)
% set plot_final==1 to plot final iterate

if c<=0 || b>=2*sqrt(c)
    error('Choose c>0 and b<2*sqrt(c)')
end

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable
x=(A/Nx)*(2*(1:Nx)-Nx); % Real spatial variable

% If no initial guess specified, use this
if u==0
    u=exp(-x.^2); 
end


% Fourier multiplier that corresponds to linear operator Lu=u''''+bu''+cu
m = pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;


% Determine alpha for which P(alpha*u)=0 and rescale u
N=sum((abs(u).^(p+1)))*(A/Nx);
dJ=real(ifft(m.*fft(u)));
J=(1/2)*sum(dJ.*u)*(A/Nx);
alpha=(2*J/N)^(1/(p-1));
u=alpha*u;


% Initialize step counter and errors

E=[];
EDiff=[];

Reserr=tol+1;


% Main loop. Iterate until H^2 norm of gradient of S is smaller than tolerance

while (Reserr>tol)

    Oldu=u;


    % Compute f(u)
    f=abs(u).^(p-1).*u; 

    % Compute gradient of S at u
    gradS=real(u-ifft(fft(f)./m)); 

    % Compute gradient of P at u
    gradP=real(2*u-(p+1)*ifft(fft(f)./m));

    % Compute the projection of gradient of S at u onto the tangent space of the
    % Nehari manifold at u by subtracting component in direction of gradient of P 

    proj = gradS-real(sum(m.*fft(gradS).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;

   

    % Gradient descent step
    u=u-delta*proj;


    % Determine alpha for which P(alpha*u)=0 and rescale u
    N=sum((abs(u).^(p+1)))*(A/Nx);
    dJ=real(ifft(m.*fft(u)));
    J=(1/2)*sum(dJ.*u)*(A/Nx);
    alpha=(2*J/N)^(1/(p-1));
    u=alpha*u;

    % Compute residual error, the H^2 norm of the gradient of S, and store in error vector
    LgradS=real(ifft(fft(gradS).*m));
    Reserr=sqrt(sum(LgradS.*gradS)*(A/Nx))
    E=[E;Reserr];

    % Compute H^2 error between iterates and store in error vector
    diff=u-Oldu;
    Ldiff=real(ifft(fft(diff).*m));
    Differr=sqrt(sum(Ldiff.*diff)*(A/Nx));
    EDiff=[EDiff; Differr];

end


% Plot final iterate
if plot_final==1
    figure(1)
    hold off
    plot(x,u,'LineWidth',1)
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
    ylabel({'$u$'},'interpreter','latex','FontSize',12)
end

