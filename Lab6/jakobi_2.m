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
n=35;

h=(xb-xa)/(n+1);
k=(yd-yc)/(xb-xa)*(n+1)-1;
x=[xa:h:xb];
y=[yc:h:yd];

tol=1e-4;
error = 1; 
licznik=0; 

%tworzenie macierzy
U1(1:n+2) = u1(x);
U2(1:k+2) = u2(y(1:(k+2)));
U3(1:n+2) = u3(x);
U4(1:k+2) = u4(y(1:(k+2)));

U(1,:) = U1(1:n+2);
U(k+2,:) = U3(1:n+2);
U(:,1) = U4(1:k+2);
U(:,n+2) = U2(1:k+2);

Uk=U;

while error>tol
    licznik = licznik+1;
    for i=2:k+1
        for j=2:n+1
            Uk(i,j) =0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1))-0.25*h^2*F(x(j),y(i));
        end
    end
    
    error = max(max(abs(Uk-U)));

    U=Uk;
end

%wykresy
[X,Y] = meshgrid(x,y);
subplot(1,2,1)
surf(X,Y,U)
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X,Y,(G(X,Y)))
title('Metoda Analityczna')

blad = max(max(abs(U-G(X,Y))))
licznik
toc
