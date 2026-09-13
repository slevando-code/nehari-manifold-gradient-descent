function ratioplot_p_exact(p,b,A,Nx,tol,u,delta)

% Script to plot error ratios for four different methods of approximating
% solutions of u''''+bu''+cu=|u|^(p-1)*u
% Algorithm 1 = Gradient descent (gradsolver_p_exact)
% Algorithm 2 = Reimannian gradient descent (rgradsolver_p_exact)
% Algorithm 3 = Inexact Newton (rnewtsolver_p_exact)
% Sample usage: ratioplot_p_exact(3,-2.5,100,4000,1e-10,0,1)

if b>=0
    error('Choose b<0.')
end

newE=gradsolver_p_exact(p,b,A,Nx,tol,u,delta,0,0);
figure(3)
hold off
plot(newE(2:end)./newE(1:end-1),'k^-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
hold on

newE=rgradsolver_p_exact(p,b,A,Nx,tol,u,delta,0,0);
figure(3)
hold on
plot(newE(2:end)./newE(1:end-1),'ko-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');

newE=rnewtsolver_p_exact(p,b,A,Nx,tol,u,delta,1,0,0);
figure(3)
hold on
plot(newE(2:end)./newE(1:end-1),'kv-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');

newE=rnewtsolver_p_exact(p,b,A,Nx,tol,u,delta,2,0,0);
figure(3)
hold on
plot(newE(2:end)./newE(1:end-1),'ks-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');

xlabel({'Iteration'},'interpreter','latex','FontSize',12)
ylabel({'$\frac{\|u_{k+1}-\varphi\|}{\|u_{k}-\varphi\|}$'},'interpreter','latex','FontSize',16)
legend('Algorithm 1','Algorithm 2','Algorithm 3 ($m=1$)','Algorithm 3 ($m=2$)','interpreter','latex','FontSize',12,'Location','best')
set(gca, 'TickLabelInterpreter', 'latex')
