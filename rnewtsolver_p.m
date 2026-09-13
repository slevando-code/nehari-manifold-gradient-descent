function E=rnewtsolver_p(p,b,c,A,Nx,tol,u0,delta,n,plot_last)

% Inexact Newton method for approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Sample usage: E=rnewtsolver_p(3,1.8,1,50,2000,1e-10,0,1,2,1)
% Spatial interval is [-A,A], with Nx subintervals. 
% delta = step size in gradient descent step
% u0 = initial guess
% tol = error tolerance 
% n = # of terms in Neumann series approximation of inverse of Hessian of action functional S
% Output E is vector of H^2 errors between iterates and exact solution

k=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable

x=(A/Nx)*(2*(1:Nx)-Nx); % Real spatial variable

% If no initial guess specified, use this
u=u0;
if u==0
    u=exp(-x.^2); 
end

m = pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;


% Determine alpha for which P(alpha*u)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=(2*I/N2)^(1/(p-1));
u=alpha*u;

Reserr=tol+1;
E=[];

while (Reserr>tol)

Oldu=u;

f=abs(u).^(p-1).*u;
df=p*abs(u).^(p-1);

%Compute gradient of S at =u
gradS=real(u-ifft(fft(f)./m));

%Compute gradient of P at u

gradP=real(2*u-ifft(fft(f+u.*df)./m));

%Project onto tangent space of N 

w0 = gradS-real(sum(m.*fft(gradS).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
w=w0;

% Approximate inverse of Hessian of S at u by truncated Neumann series

for i=1:n
   J=ifft(fft(df.*w)./m);
   J=J-real(sum(m.*fft(J).*conj(fft(gradP)))/sum(m.*fft(gradP).*conj(fft(gradP))))*gradP;
   w=w0+(1-delta)*w+delta*J;
end


u=u-delta*w;


% Determine alpha for which P(alpha*Phi)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
dI=real(ifft(m.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=(2*I/N2)^(1/(p-1));

u=alpha*u;


LgradS=real(ifft(m.*fft(gradS)));

Reserr=sqrt(sum(LgradS.*gradS)*(A/Nx));

diff=u-Oldu;

Ldiff=real(ifft(m.*fft(diff)));

differr=sqrt(sum(Ldiff.*diff)*(A/Nx));
E=[E;differr];

end


if plot_last==1
    figure(1)
    hold off
    plot(x,u)
end
