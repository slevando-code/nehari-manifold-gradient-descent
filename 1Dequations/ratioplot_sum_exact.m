function ratioplot_sum_exact(p,b,A,Nx,tol,u,delta)

% Plot ratios of H^2 norm of differences between successive iterates for
% each algorithm approximating solutions of
% u''''+bu''+cu=|u|^(p-1)*u+|u|^(q-1)*u where q=2*p-1
% Algorithm 1.1 = Gradient descent (gradsolver_sum)
% Algorithm 1.2 = Reimannian gradient descent (rgradsolver_sum)
% Algorithm 1.3 = Inexact Newton (rnewtsolver_sum)
% Sample usage: ratioplot_sum_exact(3,-4,100,4000,1e-10,0,1)

newE=gradsolver_sum_exact(p,b,A,Nx,tol,u,delta,0,0);
figure(7)
hold off
plot(newE(2:end)./newE(1:end-1),'k^-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
hold on

newE=rgradsolver_sum_exact(p,b,A,Nx,tol,u,delta,0,0);

figure(7)
hold on
plot(newE(2:end)./newE(1:end-1),'ko-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
newE=rnewtsolver_sum_exact(p,b,A,Nx,tol,u,delta,1,0,0);
figure(7)
hold on
plot(newE(2:end)./newE(1:end-1),'kv-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
newE=rnewtsolver_sum_exact(p,b,A,Nx,tol,u,delta,2,0,0);
figure(7)
hold on
plot(newE(2:end)./newE(1:end-1),'ks-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');


xlabel({'Iteration'},'interpreter','latex','FontSize',12)
ylabel({'$\frac{\|u_{k+1}-\varphi\|}{\|u_{k}-\varphi\|}$'},'interpreter','latex','FontSize',16)
legend('Algorithm 1.1','Algorithm 1.2','Algorithm 1.3 ($m=1$)','Algorithm 1.3 ($m=2$)','interpreter','latex','FontSize',12,'Location','best')
set(gca, 'TickLabelInterpreter', 'latex')
