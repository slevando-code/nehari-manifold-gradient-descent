function V =mudeltaplot_sum(p,q,a,b,c,A,B,Nx,Ny,tol,u,Delta)

% Plots asymptotic error constant mu in terms of step size delta
% Algorithm 1.1 = Gradient descent (gradsolver_sum)
% Algorithm 1.2 = Reimannian gradient descent (rgradsolver_sum)
% Algorithm 1.3 = Inexact Newton (rnewtsolver_sum)
% Sample usage:
% V=mudeltaplot_sum(3,6,1,-1.7,0,20,20,500,500,1e-10,0,.4:.1:1.8)
% V=mudeltaplot_sum(3,6,1,0,1.35,20,20,500,500,1e-10,0,.4:.1:1.8)
% V=mudeltaplot_sum(3,6,0,1,0,10,10,500,500,1e-10,0,.4:.1:1.8)

Mu1=[];
Mu2=[];
Mu3=[];
Mu4=[];

n=1;

for delta=Delta

delta
newE=gradsolver_sum(p,q,a,b,c,A,B,Nx,Ny,tol,u,delta,0);
mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
Mu1=[Mu1 mu];

newE=rgradsolver_sum(p,q,a,b,c,A,B,Nx,Ny,tol,u,delta,0);
mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
Mu2=[Mu2 mu];

newE=rnewtsolver_sum(p,q,a,b,c,A,B,Nx,Ny,tol,u,delta,1,0);
mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
Mu3=[Mu3 mu];

newE=rnewtsolver_sum(p,q,a,b,c,A,B,Nx,Ny,tol,u,delta,2,0);
mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
Mu4=[Mu4 mu];

end

V=[Delta' Mu1' Mu2' Mu3' Mu4'];

width=1;

figure(6)
hold off
plot(Delta,Mu1,'k--','LineWidth',width)
hold on
plot(Delta,Mu2,'k-','LineWidth',width)
plot(Delta,Mu3,'k:','LineWidth',width)
plot(Delta,Mu4,'k-.','LineWidth',width)

xlabel({'$\delta$'},'interpreter','latex','FontSize',12)
ylabel({'$\mu$'},'interpreter','latex','FontSize',12)
legend('Algorithm 1.1','Algorithm 1.2','Algorithm 1.3 ($m=1$)','Algorithm 1.3 ($m=2$)','interpreter','latex','FontSize',12,'Location','best')
set(gca, 'TickLabelInterpreter', 'latex')



