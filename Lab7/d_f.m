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
m=40;
D=1;
deltax=(xb-xa)/(m-1);
x=[xa:deltax:xb];         %przedzia� przestrzenny
deltat=(deltax^2)/(10*D); 
n_end=floor(yd/deltat)+1;
t=[0:deltat:1];           %przedzia� czasowy

%macierz
psi=zeros(n_end,length(x)); %utworzenie pustej macierzy

%dodanie warunku pocz�tkowego
psi(1,:) = u3(x);
psi(:,1) = u1(t);
psi(:,m) = u2(t);

psi(2,:) = psi(1,:);
F = diag(eye(2*m-2));
A1 = eye(m-2) .* -(D/deltax^2-1/(2*deltat));
A2 = (eye(m-2,m) + diag(diag(eye(m-2)),2)(1:m-2,:)) .* D/deltax^2;
A = [A1, A2] ./ (1/(2*deltat) + D/deltax^2);
for n=3:n_end
  F(1:m-2) = psi(n-2,2:m-1);
  F(m-1:length(F)) = psi(n-1,:);
  psi(n,2:m-1) = A * F;
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