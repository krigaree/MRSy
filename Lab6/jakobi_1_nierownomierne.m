%metoda Jacobiego
clc
clear all
tic

%funkcja
F = @(x,y) 0;

%rozwi�zanie analityczne
G = @(x,y) log(x.^2+y.^2);

%przedzia� omega
xa=1;
xb=2;
yc=0;
yd=1;

%warunki brzegowe
war1 = @(x) 2*log(x);
war2 = @(y) log(y.^2+4);
war3 = @(x) log(x.^2+1);
war4 = @(y) log(y.^2+1);

%siatka
n=25;
m=25;

h=(xb-xa)/(n+1);
k=(yd-yc)/(m+1);
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

for i=1:m+2
    for j=1:n+2
        g(i,j) = G(x(j),y(i));
    end
end

while error>tol
    for i=2:m+1
        for j=2:n+1
            M(i-1,j-1) = (1/h^2*(M0(i+1,j)+M0(i-1,j))+1/k^2*(M0(i,j+1)+M0(i,j-1)));
            M(i-1,j-1) = M(i-1,j-1) - F(x(j),y(i));
            M(i-1,j-1) = M(i-1,j-1)/(2*(1/h^2 + 1/k^2));
        end
    end
    
    error=max(max(abs(M0(2:m+1,2:n+1)-M)));
    M0(2:m+1,2:n+1)=M;
    licznik = licznik+1;
    
end

toc
%wykresy
[X,Y] = meshgrid(x,y);
subplot(1,2,1)
surf(X,Y,M0)
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X,Y,(G(X,Y)))
title('Metoda Analityczna')
Error=max(max(abs(M0-G(X,Y))))
licznik

