function mudeltaplot_p_exact(p,b,A,Nx,tol,u,Delta)

% Script for plotting asymptotic error constant mu in terms of step size
% delta for four methods of approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Algorithm 1 = Gradient descent (gradsolver_p_exact)
% Algorithm 2 = Reimannian gradient descent (rgradsolver_p_exact)
% Algorithm 3 = Inexact Newton (rnewtsolver_p_exact) 

Mu1=[];
Mu2=[];
Mu3=[];
Mu4=[];

for delta=Delta

E=gradsolver_p_exact(p,b,A,Nx,tol,u,delta,0,0);
%figure(7)
%hold off
%plot(newE(2:end)./newE(1:end-1),'k^-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
%hold on

r=E(end-3:end)./E(end-4:end-1);
mu=mean(r);
Mu1=[Mu1 mu];

E=rgradsolver_p_exact(p,b,A,Nx,tol,u,delta,0,0);
r=E(end-3:end)./E(end-4:end-1);
mu=mean(r);
Mu2=[Mu2 mu];

E=rnewtsolver_p_exact(p,b,A,Nx,tol,u,delta,1,0,0);
r=E(end-3:end)./E(end-4:end-1);
mu=mean(r);
Mu3=[Mu3 mu];

E=rnewtsolver_p_exact(p,b,A,Nx,tol,u,delta,2,0,0);
r=E(end-3:end)./E(end-4:end-1);
mu=mean(r);
Mu4=[Mu4 mu];

end

E=petsolver_p_exact(p,b,A,Nx,tol,u,0,0);
r=E(end-3:end)./E(end-4:end-1);
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
legend('Algorithm 1','Algorithm 2','Algorithm 3 ($m=1$)','Algorithm 3 ($m=2$)','Petviashvili','interpreter','latex','FontSize',12,'Location','best')
set(gca, 'TickLabelInterpreter', 'latex')
