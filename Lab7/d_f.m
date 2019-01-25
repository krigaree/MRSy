clc,clear all
tic
%rozwi�zanie analityczne
G = @(x,t) sin(pi.*x./2).*exp(-(pi.^2).*t./4);

%przedzia� omega
xa=0;
xb=2;
yc=0;
yd=1;

%warunki brzegowe
u1 = @(x) 0;
u2 = @(x) 0;
u3 = @(x,t) sin(pi*x/2);

licznik=0;
%siatka
m=20;
D=1;
deltax=(xb-xa)/(m-1);
x=[xa:deltax:xb];         %przedzia� przestrzenny
deltat=(deltax^2)/(20*D); %dzielimy od razu przez 10, aby warto�� nie by�a blisko deltat graniczne
n_end=floor(yd/deltat)+1;
t=[0:deltat:1];           %przedzia� czasowy

%macierz
psi=zeros(n_end,length(x)); %utworzenie pustej macierzy

%dodanie warunku pocz�tkowego
psi(1,:) = u3(x);
psi(:,1) = u1(t);
psi(:,m) = u2(t);

%A = D/deltax^2 *diag(eye(m-2));
%B = diag(A) + D/deltax^2*diag(diag(eye(m-3)),-1) + -1*diag(diag(eye(m-3)),1);
psi(2,:) = psi(1,:);
for n=3:n_end
  for i=2:m-1
    psi(n,i) = (D/deltax^2) * (psi(n-1,i+1) - psi(n-2,i) + psi(n-1,i-1)) + psi(n-2,i)/(2*deltat);
    psi(n,i) = psi(n,i) / (1/(2*deltat) + D/deltax^2);
  end
  licznik = licznik+1;
end

[X,T] = meshgrid(x,t);
subplot(1,2,1)
surf(X,T,psi)
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X,T,(G(X,T)))
title('Metoda Analityczna')
Error=max(max(abs(psi-G(X,T))));
licznik
toc