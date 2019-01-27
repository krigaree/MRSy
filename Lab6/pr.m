%metoda Jacobiego
clc
clear all
%tic

%funkcja
F = @(x,y) 0;

%rozwi?zanie analityczne
G = @(x,y) log(x.^2+y.^2);

%przedzia? omega
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
n=30;
m=30;
h=(xb-xa)/(n+1);
k=(yd-yc)/(m+1);
x=[xa:h:xb];
y=[yc:k:yd];

tol=1e-4;
error = 1; 
licznik=0; 

%tworzenie macierzy
U1(1:n+2) = u1(x);
U2(1:m+2) = u2(y(1:(m+2)));
U3(1:n+2) = u3(x);
U4(1:m+2) = u4(y(1:(m+2)));

U(1,:) = U1(1:n+2);
U(m+2,:) = U3(1:n+2);
U(:,1) = U4(1:m+2);
U(:,n+2) = U2(1:m+2);

Uk=U;
R = createR(n+1);
lR = length(R); 
while error>tol
      r = R(mod(licznik,lR)+1);
      %step n -> n+1/2
      step = true;
      iter = length(y);
      for i=2:iter-1
        Uk(i, 2:length(x)-1) = doStep(i, r, F, x, y, U, step);
      end
      U = Uk;
      %step n+1/2 -> n+1
      step = false; 
      iter = length(x);
      for i=2:iter-1
        Uk(2:length(y)-1, i) = doStep(i, r, F, x, y, U, step);
      end
      licznik = licznik+1;
      error = max(max(abs(Uk-U)));
      Error(licznik) = error;

      U=Uk;
end

%wykresy
n = 100;
m = 100;
h2=(xb-xa)/(n+1);
k2=(yd-yc)/(m+1);
x2=[xa:h2:xb];
y2=[yc:k2:yd];
[X,Y] = meshgrid(x,y);
[X2,Y2] = meshgrid(x2,y2);
subplot(1,2,1)
surf(X,Y,U)
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X2,Y2,(G(X2,Y2)))
title('Metoda Analityczna')

blad = max(max(abs(U-G(X,Y))))
licznik
%toc