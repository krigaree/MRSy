%metoda Gaussa - Seidela
clc
clear all
tic

%funkcja
F = @(x,y) -cos(x+y)-cos(x-y);

%rozwi�zanie analityczne
G = @(x,y) cos(x).*cos(y);

%przedzia� omega
xa=0;
xb=pi;
yc=0;
yd=pi/2;

%warunki brzegowe
war1 = @(x) cos(x);
war2 = @(y) -cos(y);
war3 = @(x) 0;
war4 = @(y) cos(y);

%siatka
n=50;
m=(yd-yc)/(xb-xa)*(n+1)-1;

h=(xb-xa)/(n+1);
x=linspace(xa,xb,n+2);
y=linspace(yc,yd,m+2);

tol=1e-4; 
error = 10; 
licznik=0; 

%tworzenie macierzy
M0=zeros(m,n)+1;
M1(1:n+2) = war1(x);
M2(1:m) = war2(y(2:length(y)-1));
M3(1:n+2) = war3(x);
M4(1:m) = war4(y(2:length(y)-1));

M0=[M4', M0, M2'];
M0=[M1;M0;M3];
M=M0;

for i=1:m+2
    for j=1:n+2
        g(i,j) = G(x(j),y(i));
    end
end

while error>tol
    for i=2:m+1
        for j=2:n+1
            M(i,j) =0.25*(M(i+1,j)+M(i-1,j)+M(i,j+1)+M(i,j-1))-0.25*h^2*F(x(j),y(i));
        end
    end

    error=max(max(abs(M0-M)));
    M0=M;
    licznik = licznik+1;
    
end

%wykresy
[X,Y] = meshgrid(x,y);
subplot(1,2,1)
surf(X,Y,M0)
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X,Y,(G(X,Y)))
title('Metoda Analityczna')

licznik
toc
