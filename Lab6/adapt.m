%metoda Peacemanna-Rachforda
tic
clc
clear all

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
i=1;
n=i;
m=n;
Ut = zeros(m+2,m+2);
x=0;
y=0;
U=0;
tol=1e-4;
error = 1; 
licznik=0;
%siatka
krok = 3;
while krok < 100 %krok = dlugosc siatki krok x krok xD
  licznik = 0;
  krok = krok*2-1;
  n=krok-2;
  m=n;
  h=(xb-xa)/(n+1);
  k=(yd-yc)/(m+1);
  x=[xa:h:xb];
  y=[yc:k:yd];

  %tworzenie macierzy
  U1 = u1(x);
  U2 = u2(y);
  U3 = u3(x);
  U4 = u4(y);

  U = zeros(size(Ut)*2-1);
  U(1,:) = U1;
  U(m+2,:) = U3;
  U(:,1) = U4;
  U(:,n+2) = U2;
  if i >3
    U(3:2:length(U(:,1))-2, 3:2:length(U(1,:))-2) = Ut(2:length(Ut(:,1))-1, 2:length(Ut(1,:))-1);
  end
  Ut = U;
  for i=2:length(U(:,1))-1
    for j=2:length(U(1,:))
      if (U(i,j) == 0)
        if nnz(Ut(i-1:i+1,j-1:j+1)) > 0
          U(i,j) = sum(sum(Ut(i-1:i+1,j-1:j+1))) / nnz(Ut(i-1:i+1,j-1:j+1));
        end
      end
    end  
  end
  U(1,:) = U1;
  U(m+2,:) = U3;
  U(:,1) = U4;
  U(:,n+2) = U2;

  Uk=U;
  mR = max(n,m);
  R = createR(mR+1);
  lR = length(R); 
  while licznik < 5
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
        %Error(licznik) = error;

        U=Uk;
  end
  licznik = licznik +1;
  Ut = U;
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