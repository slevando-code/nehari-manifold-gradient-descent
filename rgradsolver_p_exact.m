function E=rgradsolver_p_exact(p,b,A,Nx,tol,u,delta,plot_last,plot_error)

% Riemannian Nehari gradient descent method for approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Sample usage: E=rgradsolver_p_exact(3,-2.5,50,2000,1e-10,0,1.2,1,1)
% Spatial interval is [-A,A], with Nx subintervals. 
% delta = step size in gradient descent step
% u = initial guess
% tol = error tolerance 
% b must be less than 0
% Output E is vector of H^2 errors between iterates and exact solution

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



% Determine alpha for which P(alpha*Phi)=0 and rescale u
N=sum((abs(u).^(p+1)))*(A/Nx);
dJ=real(ifft(m.*fft(u)));
J=(1/2)*sum(dJ.*u)*(A/Nx);
alpha=(2*J/N)^(1/(p-1));
u=alpha*u;

% Initialize counter and errors

H2err=tol+1;
E=[];

while (H2err>tol)

Oldu=u;


f=abs(u).^(p-1).*u;

gradS=real(u-ifft(fft(f)./m));

gradP=real(2*u-(p+1)*ifft(fft(f)./m));


proj = gradS-real(sum(m.*fft(gradS).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;

u=u-delta*proj;


% Determine alpha for which P(alpha*Phi)=0
N=sum((abs(u).^(p+1)))*(A/Nx);
dJ=real(ifft(m.*fft(u)));
J=(1/2)*sum(dJ.*u)*(A/Nx);
alpha=(2*J/N)^(1/(p-1));
u=alpha*u;



Exactdiff=u-exact;
LExactdiff=real(ifft(fft(Exactdiff).*m));
H2err=sqrt(sum(LExactdiff.*Exactdiff)*(A/Nx));

E=[E; H2err];

end



% Plot exact solution together with final iterate
if plot_last==1
    figure(1)
    hold off
    plot(x,u,'b')
    hold on
    plot(x,exact,'r')
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
    ylabel({'$\varphi(x)$'},'interpreter','latex','FontSize',12)
end


% Plot ratio of successive H^2 errors between iterate and exact solution.
%figure(2)
%hold off
%plot(E(2:end)./E(1:end-1))
%hold on
%plot(ERes(2:end)./ERes(1:end-1))
%plot(EDiff(2:end)./EDiff(1:end-1))
%xlabel({'Iteration'},'interpreter','latex','FontSize',12)
%legend({'$\frac{\|u_{k+1}-\varphi\|}{\|u_k-\varphi\|}$','$\frac{\|\nabla S(u_{k+1})\|}{\|\nabla S(u_k)\|}$','$\frac{\|u_{k+2}-u_{k+1}\|}{\|u_{k+1}-u_k\|}$'},'interpreter','LaTex','FontSize',12)


% Plot log of exact errors, residual errors and difference between iterates
%figure(3)
%hold off
%plot(log(E))
%hold on
%plot(log(ERes))
%plot(log(EDiff))
%legend({'$\log\|u_k-\varphi\|$','$\log\|\nabla S(u_k)\|$','$\log\|u_{k+1}-u_k\|$'},'interpreter','LaTex','FontSize',12)

%xlabel({'Iteration'},'interpreter','latex','FontSize',12)


% Plot just ratio of successive H^2 errors between iterate and exact solution.
%figure(4)
%hold off
%plot(E(2:end)./E(1:end-1),'k-o','MarkerSize',3,'MarkerFaceColor','k')
%plot(E(2:end)./E(1:end-1),'k-o','MarkerSize',3,'MarkerFaceColor','k')
%hold on
%xlabel({'Iteration'},'FontSize',12)
%xlabel({'Iteration'},'interpreter','latex','FontSize',16)
%ylabel({'$\frac{\|u_{k+1}-\varphi\|}{\|u_k-\varphi\|}$'},'interpreter','latex','FontSize',16)
%s=strcat('$\mu\approx$',num2str(sigma));
%legend({'$\frac{\|u_{k+1}-\varphi\|}{\|u_k-\varphi\|}$'},'interpreter','LaTex','FontSize',20)
%title(s,'interpreter','LaTex','FontSize',16)

% Plot exact error of last iterate
if plot_error==1
    figure(5)
    plot(x,Exactdiff,'k')
    xlabel({'$x$'},'interpreter','latex','FontSize',12)
end

