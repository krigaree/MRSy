%metoda Jacobiego
clc
clear all
tic

%funkcja
F = @(x,y) 0;

%rozwi¹zanie analityczne
G = @(x,y) log(x.^2+y.^2);

%przedzia³ omega
xa=1;
xb=2;
yc=0;
yd=1;

%warunki brzegowe
u1 = @(x) 2*log(x);
u2 = @(y) log(y.^2+4);
u3 = @(x) log(x.^2+1);
u4 = @(y) log(y.^2+1);

%siatka
n=50;

h=(xb-xa)/(n+1);
k=(yd-yc)/(n+1);
x=[xa:h:xb];
y=[yc:k:yd];

tol=1e-4; 
error = 1; 
licznik=0; 

%tworzenie macierzy
for i=1:n+2
    U(i,1)=2*log(xa+(i-1)*h);
    U(1,i)=log(((yc+(i-1)*k)^2)+1);
end

for i=2:n+2
    U(i,n+2)=log(((xa+(i-1)*h)^2)+1);
    U(n+2,i)=log(((yc+(i-1)*k)^2)+4);
end

Uk=U;

while error>tol
    licznik = licznik+1;
    for i=2:n+1
        for j=2:n+1
            Uk(i,j) =0.25*(U(i+1,j)+U(i-1,j)+U(i,j+1)+U(i,j-1))-0.25*h^2*F(x(i),y(j));
        end
    end
    
    error = max(max(abs(Uk-U)));

    U=Uk;
end

%wykresy
[X,Y] = meshgrid(x,y);
subplot(1,2,1)
surf(X,Y,U')
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X,Y,(G(X,Y)))
title('Metoda Analityczna')

blad = max(max(abs(U'-G(X,Y))))
licznik
toc
