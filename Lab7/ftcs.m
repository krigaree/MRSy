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
m=10;
D=1;
deltax=(xb-xa)/(m-1);
x=[xa:deltax:xb];         %przedzia� przestrzenny
deltat=(deltax^2)/(20*D); %dzielimy od razu przez 10, aby warto�� nie by�a blisko deltat graniczne
n_end=floor(yd/deltat)+1;
t=[0:deltat:1];           %przedzia� czasowy

%macierz
psi=zeros(n_end,length(x)); %utworzenie pustej macierzy

for i=2:m-1                %dodanie warunku pocz�tkowego
  psi(1,i)=u3(x(i));
end
                        
for n=2:n_end
  for i=2:m-1
    psi(n,i)=((D*deltat)/deltax^2)*(psi(n-1,i+1)-2*psi(n-1,i)+psi(n-1,i-1))+psi(n-1,i);
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
Error=max(max(abs(psi-G(X,T))))
licznik;
toc