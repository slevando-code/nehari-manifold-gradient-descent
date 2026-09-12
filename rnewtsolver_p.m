function EDiff=rnewtsolver_p(p,b,c,A,Nx,tol,u0,delta,m,plot_last)

% Inexact Newton method for approximating solutions of u''''+bu''+cu=|u|^(p-1)*u
% Sample usage rgradsolver_p(3,-2,1,40,500,1e-10,0,1)
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

v= pi^4*k.^4/A^4-b*pi^2*k.^2/A^2+c;


% Determine alpha for which P(alpha*u)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
dI=real(ifft(v.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=(2*I/N2)^(1/(p-1));
u=alpha*u;

Reserr=tol+1;
E=[];
EDiff=[];

while (Reserr>tol)

Oldu=u;

f=abs(u).^(p-1).*u;
df=p*abs(u).^(p-1);

%Compute gradient of S at =u
gradS=real(u-ifft(fft(f)./v));

%Compute gradient of P at u

gradP=real(2*u-ifft(fft(f+u.*df)./v));

%Project onto tangent space of N 

w0 = gradS-real(sum(v.*fft(gradS).*conj(fft(gradP)))/sum(v.*fft(gradP).*conj(fft(gradP))))*gradP;
w=w0;

% Approximate inverse of Hessian of S at u by truncated Neumann series

for i=1:m
   J=ifft(fft(df.*w)./v);
   J=J-real(sum(v.*fft(J).*conj(fft(gradP)))/sum(v.*fft(gradP).*conj(fft(gradP))))*gradP;
   w=w0+(1-delta)*w+delta*J;
end


u=u-delta*w;


% Determine alpha for which P(alpha*Phi)=0
N2=sum((abs(u).^(p+1)))*(A/Nx);
dI=real(ifft(v.*fft(u)));
I=(1/2)*sum(dI.*u)*(A/Nx);
alpha=(2*I/N2)^(1/(p-1));

u=alpha*u;


LgradS=real(ifft(v.*fft(gradS)));

Reserr=sqrt(sum(LgradS.*gradS)*(A/Nx));
E=[E;Reserr];


diff=u-Oldu;

Ldiff=real(ifft(v.*fft(diff)));

differr=sqrt(sum(Ldiff.*diff)*(A/Nx));
EDiff=[EDiff;differr];

end


if plot_last==1
    figure(1)
    hold off
    plot(x,u)
end
