function V=mudeltaplot_p_exact(p,b,A,Nx,tol,u,Delta)

% Script for plotting asymptotic error constant mu in terms of step size
% delta for four methods of approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Algorithm 1.1 = Gradient descent (gradsolver_p_exact)
% Algorithm 1.2 = Reimannian gradient descent (rgradsolver_p_exact)
% Algorithm 1.3 = Inexact Newton (rnewtsolver_p_exact)
% Sample usage: V=mudeltaplot_p_exact(3,-2.5,100,2000,1e-10,0,.4:.01:1.8)

Mu1=[];
Mu2=[];
Mu3=[];
Mu4=[];

n=1;

for delta=Delta

    delta

    E=gradsolver_p_exact(p,b,A,Nx,tol,u,delta,0,0);
    r=E(end-n:end)./E(end-n-1:end-1);
    mu=mean(r);
    Mu1=[Mu1 mu];

    E=rgradsolver_p_exact(p,b,A,Nx,tol,u,delta,0,0);
    r=E(end-n:end)./E(end-n-1:end-1);
    mu=mean(r);
    Mu2=[Mu2 mu];

    E=rnewtsolver_p_exact(p,b,A,Nx,tol,u,delta,1,0,0);
    r=E(end-n:end)./E(end-n-1:end-1);
    mu=mean(r);
    Mu3=[Mu3 mu];

    E=rnewtsolver_p_exact(p,b,A,Nx,tol,u,delta,2,0,0);
    r=E(end-n:end)./E(end-n-1:end-1);
    mu=mean(r);
    Mu4=[Mu4 mu];

end

V=[Delta' Mu1' Mu2' Mu3' Mu4'];

E=petsolver_p_exact(p,b,A,Nx,tol,u,0,0);
r=E(end-n:end)./E(end-n-1:end-1);
petmu=mean(r);

w=1;

figure(6)
hold off
plot(Delta,Mu1,'k--','LineWidth',w)
hold on
plot(Delta,Mu2,'k-','LineWidth',w)
plot(Delta,Mu3,'k:','LineWidth',w)
plot(Delta,Mu4,'k-.','LineWidth',w)

plot(1,petmu,'k.','MarkerSize',12)


xlabel({'$\delta$'},'interpreter','latex','FontSize',12)
ylabel({'$\mu$'},'interpreter','latex','FontSize',12)
legend('Algorithm 1.1','Algorithm 1.2','Algorithm 1.3 ($m=1$)','Algorithm 1.3 ($m=2$)','Petviashvili','interpreter','latex','FontSize',12,'Location','best')
set(gca, 'TickLabelInterpreter', 'latex')
