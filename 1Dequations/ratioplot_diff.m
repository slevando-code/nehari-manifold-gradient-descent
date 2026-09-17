function ratioplot_diff(p,q,b,c,A,Nx,tol,u,delta)

% Plot ratios of H^2 norm of differences between successive iterates for
% each algorithm approximating solutions of u''''+bu''+cu=|u|^(q-1)*u-|u|^(p-1)*u
% Algorithm 1.1 = Gradient descent (gradsolver_sum)
% Algorithm 1.2 = Reimannian gradient descent (rgradsolver_sum)
% Algorithm 1.3 = Inexact Newton (rnewtsolver_sum)
% Sample usage: ratioplot_diff(3,6,-3,1,100,4000,1e-10,0,1)

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2];
x=(A/Nx)*(2*(1:Nx)-Nx);

if u==0
    u=exp(-(x.^2));
end

newE=gradsolver_diff(p,q,b,c,A,Nx,tol,u,delta,0);
figure(7)
hold off
plot(newE(2:end)./newE(1:end-1),'k^-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
hold on

newE=rgradsolver_diff(p,q,b,c,A,Nx,tol,u,delta,0);

figure(7)
hold on
plot(newE(2:end)./newE(1:end-1),'ko-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
newE=rnewtsolver_diff(p,q,b,c,A,Nx,tol,u,delta,1,0);
figure(7)
hold on
plot(newE(2:end)./newE(1:end-1),'kv-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');
newE=rnewtsolver_diff(p,q,b,c,A,Nx,tol,u,delta,2,0);
figure(7)
hold on
plot(newE(2:end)./newE(1:end-1),'ks-','LineWidth',.7,'MarkerSize',4,'MarkerFaceColor','k');


xlabel({'Iteration'},'interpreter','latex','FontSize',12)
ylabel({'$\frac{\|u_{k+2}-u_{k+1}\|}{\|u_{k+1}-u_k\|}$'},'interpreter','latex','FontSize',16)
legend('Algorithm 1.1','Algorithm 1.2','Algorithm 1.3 ($m=1$)','Algorithm 1.3 ($m=2$)','interpreter','latex','FontSize',12,'Location','best')
set(gca, 'TickLabelInterpreter', 'latex')
