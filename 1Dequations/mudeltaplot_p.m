function V=mudeltaplot_p(p,b,c,A,Nx,tol,u,Delta)

% Script for plotting asymptotic error constant mu in terms of step size
% delta for four methods of approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Algorithm 1.1 = Gradient descent (gradsolver_p)
% Algorithm 1.2 = Reimannian gradient descent (rgradsolver_p)
% Algorithm 1.3 = Inexact Newton (rnewtsolver_p)
% Sample usage: V=mudeltaplot_p(3,1,1,100,4000,1e-10,0,.4:.01:1.8)

Mu1=[];
Mu2=[];
Mu3=[];
Mu4=[];

n=1;

for delta=Delta

    delta

    newE=gradsolver_p(p,b,c,A,Nx,tol,u,delta,0);
    mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
    Mu1=[Mu1 mu];

    newE=rgradsolver_p(p,b,c,A,Nx,tol,u,delta,0);
    mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
    Mu2=[Mu2 mu];

    newE=rnewtsolver_p(p,b,c,A,Nx,tol,u,delta,1,0);
    mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
    Mu3=[Mu3 mu];
    
    newE=rnewtsolver_p(p,b,c,A,Nx,tol,u,delta,2,0);
    mu=mean(newE(end-n:end)./newE(end-n-1:end-1));
    Mu4=[Mu4 mu];

end

V=[Delta' Mu1' Mu2' Mu3' Mu4'];

w=1;

figure(6)
hold off
plot(Delta,Mu1,'k--','LineWidth',w)
hold on
plot(Delta,Mu2,'k-','LineWidth',w)
plot(Delta,Mu3,'k:','LineWidth',w)
plot(Delta,Mu4,'k-.','LineWidth',w)

xlabel({'$\delta$'},'interpreter','latex','FontSize',12)
ylabel({'$\mu$'},'interpreter','latex','FontSize',12)
legend('Algorithm 1.1','Algorithm 1.2','Algorithm 1.3 ($m=1$)','Algorithm 1.3 ($m=2$)','interpreter','latex','FontSize',12,'Location','best')
set(gca, 'TickLabelInterpreter', 'latex')
