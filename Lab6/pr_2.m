%metoda Peacemanna-Rachforda
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
u1 = @(x) cos(x);
u2 = @(y) -cos(y);
u3 = @(x) 0;
u4 = @(y) cos(y);

%siatka
n=50;
m=n;
h=(xb-xa)/(n+1);
k=(yd-yc)/(m+1);
x=[xa:h:xb];
y=[yc:k:yd];

tol=1e-7;
error = 1; 
licznik=0; 

%tworzenie macierzy
U1 = u1(x);
U2 = u2(y);
U3 = u3(x);
U4 = u4(y);

U = zeros(m+2,n+2);
U(1,:) = U1;
U(m+2,:) = U3;
U(:,1) = U4;
U(:,n+2) = U2;

Uk=U;
if xb-xa > yd-yc
  mR = n;
else
  mR = m;
end
R = createR(mR+1);
lR = length(R); 
while error>tol
      r = R(mod(licznik,lR)+1);
      %step n -> n+1/2
      step = true;
      iter = length(y);
      for i=2:iter-1
        Uk(i, 2:length(x)-1) = doStep(i, r, x, y, U, step);
      end
      U = Uk;
      %step n+1/2 -> n+1
      step = false; 
      iter = length(x);
      for i=2:iter-1
        Uk(2:length(y)-1, i) = doStep(i, r, x, y, U, step);
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
surf(X,Y,(G(X,Y)))
title('Metoda Analityczna')

licznik
toc