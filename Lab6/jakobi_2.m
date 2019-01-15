%metoda Jacobiego
clc
clear all
tic

%funkcja
F = @(x,y) -cos(x+y)-cos(x-y);

%rozwi¹zanie analityczne
G = @(x,y) cos(x).*cos(y);

%przedzia³ omega
xa=0;
xb=pi;
yc=0;
yd=pi/2;

%warunki brzegowe
u1 = @(x) cos(x);
u2 = @(y) -cos(y);
u3 = @(x) 0;
u4 = @(y) cos(y);

%siatka
n=55;
m=(yd-yc)/(xb-xa)*(n+1)-1;

h=(xb-xa)/(n+1);
x=linspace(xa,xb,n+2);
y=linspace(yc,yd,m+2);

tol=1e-4; 
error = 10; 
licznik=0; 

%tworzenie macierzy
U0=zeros(m,n)+1;
U1(1:n+2) = u1(x);
U2(1:m) = u2(y(2:length(y)-1));
U3(1:n+2) = u3(x);
U4(1:m) = u4(y(2:length(y)-1));

U0=[U4', U0, U2'];
U0=[U1;U0;U3];

for i=1:m+2
    for j=1:n+2
        g(i,j) = G(x(j),y(i));
    end
end

while error>tol
    for i=2:m+1
        for j=2:n+1
            U(i-1,j-1) =0.25*(U0(i+1,j)+U0(i-1,j)+U0(i,j+1)+U0(i,j-1))-0.25*h^2*F(x(j),y(i));
        end
    end
    
    error=max(max(abs(U0(2:m+1,2:n+1)-U)));
    U0(2:m+1,2:n+1)=U;
    licznik = licznik+1;
    
end

%wykresy
[X,Y] = meshgrid(x,y);
subplot(1,2,1)
surf(X,Y,U0)
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X,Y,(G(X,Y)))
title('Metoda Analityczna')

licznik
toc
