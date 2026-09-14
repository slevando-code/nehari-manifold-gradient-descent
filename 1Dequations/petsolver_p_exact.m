function E=petsolver_p_exact(p,b,A,Nx,tol,u,plot_last,plot_error)

% Petviashvili method for approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Sample usage: E=petsolver_p_exact(3,-2.5,50,2000,1e-10,0,1,1)
% Spatial interval is [-A,A], with Nx subintervals. 
% delta = step size in gradient descent step
% u = initial guess
% tol = error tolerance 
% b must be less than 0 (for known exact solutions)
% Output E is vector of H^2 errors between iterates and exact solution

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

if ~isequal(size(x),size(u))
    u=exp(-x.^2);
end

% Fourier multiplier for the operator Lu=u''''+bu''+cu
v= pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;

gamma=p/(p-1); % optimal power of stabilizing factor 

% Initialize errors

H2err=tol+1;
E=[];

% Main Loop
while (H2err>tol)

    % Compute stabilization factor 

    f=abs(u).^(p-1).*u;
    M=real(sum(v.*fft(u).*conj(fft(u)))/sum(fft(f).*conj(fft(u))));

    % Apply Petviashvili method

    u=real(ifft((M^gamma)*fft(f)./v));

   
    % Compute H^2 error betweeen u and exact solution
    Exactdiff=u-exact;
    LExactdiff=real(ifft(fft(Exactdiff).*v));
    H2err=sqrt(sum(LExactdiff.*Exactdiff)*(A/Nx));

    % Add new H^2 error to list of errors
    E=[E; H2err];
end

% Plot exact solution together with final iterate
if plot_last==1
    figure(1)
    hold off
    plot(x,u,'b')
    hold on
    plot(x,exact,'r','LineWidth',1)
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
    ylabel({'$u$'},'interpreter','latex','FontSize',12)
end

%Plot exact error of last iterate
if plot_error==1
    figure(2)
    plot(x,Exactdiff,'k')
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
end