function E=rgradsolver_sum(p,q,a,b,c,A,B,Nx,Ny,tol,u,delta,plot_final)


% Numerically approximate solutions of
% u_xxxx+2u_xxyy+u_yyyy-b(u_xx+u_yy)+c^2u_xx+au=|u|^(p-1)*u+|u|^(q-1)*u
% Sample usage: 
% E=rgradsolver_sum(3,6,0,1,0,10,10,500,500,1e-10,0,1,1)
% E=rgradsolver_sum(3,6,1,-1.7,0,20,20,1000,1000,1e-10,0,1,1)
% E=rgradsolver_sum(3,6,1,0,1.35,20,20,1000,1000,1e-10,0,1,1)

if a<0 || c^2>=b+2*sqrt(a)
    error('Choose a>=0 and c^2<b+2*sqrt(a).')
end


kx=[(0:Nx/2) (1:Nx/2-1)-Nx/2]; % Fourier transform variable xi
x=(A/(Nx))*(2*(1:Nx)-Nx); % Real spatial variable x
ky=[(0:Nx/2) (1:Ny/2-1)-Ny/2]; % Fourier transform variable eta
y=(B/(Ny))*(2*(1:Ny)-Ny); % Real spatial variable y

[Kx,Ky]=meshgrid(kx,ky);

[X,Y]=meshgrid(x,y);

if u==0
    u=exp(-(X).^2-(Y).^2);
end


m = a*pi^4*(Kx.^2/A^2+Ky.^2/B^2).^2+b*pi^2*(Kx.^2/A^2+Ky.^2/B^2)-c^2*pi^2*Kx.^2/A^2+1;



% Determine alpha for which P(alpha*u)=0
N1=sum(sum((abs(u).^(q+1))))*(A/Nx)*(B/Ny);
N2=sum(sum((abs(u).^(p+1))))*(A/Nx)*(B/Ny);
dI=real(ifft2(m.*fft2(u)));
I=(1/2)*sum(sum(dI.*u))*(A/Nx)*(B/Ny);
alpha=PRoot_sum(p,q,2*I,N2,N1);
u=alpha*u;


E=[];
H2err=tol+1;

while (H2err>tol)

    Oldu=u;
    f=abs(u).^(p-1).*u+abs(u).^(q-1).*u;
    df=q*abs(u).^(q-1)+p*abs(u).^(p-1);

    gradS=real(u-ifft2(fft2(f)./m));

    % Compute gradient of P at u
    gradS=real(u-ifft2(fft2(f)./m));

    % Compute gradient of P at u
    gradP=real(2*u-ifft2(fft2(f+df.*u)./m));

    % Compute the projection of gradient of S at u onto the tangent space of the
    % Nehari manifold at u by subtracting component in direction of gradient of P
     proj = gradS-real(sum(sum(m.*fft2(gradS).*conj(fft2(gradP))))/sum(sum(m.*fft2(gradP).*conj(fft2(gradP)))))*gradP;

    % Gradient descent step
    u=u-delta*proj;
  
    % Determine alpha for which P(alpha*u)=0
    N1=sum(sum((abs(u).^(q+1))))*(A/Nx)*(B/Ny);
    N2=sum(sum((abs(u).^(p+1))))*(A/Nx)*(B/Ny);
    dI=real(ifft2(m.*fft2(u)));
    I=(1/2)*sum(sum(dI.*u))*(A/Nx)*(B/Ny);
    alpha=PRoot_sum(p,q,2*I,N2,N1);
    u=alpha*u;


    LgradS=real(ifft2(m.*fft2(gradS)));
    H2err=sqrt(sum(sum(LgradS.*gradS))*(A/Nx)*(B/Ny));

    diff=u-Oldu;
    Ldiff=real(ifft2(m.*fft2(diff)));
    differr=sqrt(sum(sum(Ldiff.*diff))*(A/Nx)*(B/Ny));

    E=[E;differr];

end


if plot_final==1

    figure(1)
    hold off
    s=surf(X,Y,u);


    s.EdgeColor = 'none';
    view(3)
    l1 = light;
    l1.Position = [160 400 80];
    l1.Style = 'local';
    l1.Color = [0.8 0.8 0.8];

    l2 = light;
    l2.Position = [.5 -1 .4];
    l2.Color = [0.8 0.8 0.8];
    s.FaceColor = [0.5 0.5 0.5];

    s.FaceLighting = 'gouraud';
    s.AmbientStrength = 0.3;
    s.DiffuseStrength = 0.6;
    s.BackFaceLighting = 'lit';

    s.SpecularStrength = .5;
    s.SpecularColorReflectance = .5;
    s.SpecularExponent = 5;

end

%Newton's Method for approximating the positive root of f(x)=a-b*x^(p-1)-c*x^(q-1)
function x=PRoot_sum(p,q,a,b,c)
tol=1e-15;
x=(a/b)^(1/(p-1));
e=1;

while (abs(e)>tol)
    f=a-b*x^(p-1)-c*x^(q-1);
    df=-b*(p-1)*x^(p-2)-c*(q-1)*x^(q-2);
    e=f/df;
    x=x-e;
end
