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
n=7;

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
R = createR(n+1);
lR = length(R); 
while licznik<150
    

    
    r = R(mod(licznik,lR)+1);
    %if (mod(licznik,2) == 1)
      %disp('step n')
      step = true; %step n -> n+1/2
      iter = length(y);
      for i=2:iter-1
         %start from 1 to 3
        Uk(2:iter-1, i) = doStep(i, r, F, x, y, U, step);
        Uk(2:iter-1, i);
      end
      U = Uk;
    %else %step n+1/2 -> n
      %disp('step n1/2')
      step = false;
      iter = length(x);
      for i=2:iter-1
        %r = R(mod(i+1,lR)+1); %start from 1 to 3
        Uk(i, 2:iter-1) = doStep(i, r, F, x, y, U, step);
        Uk(i, 2:iter-1);
      end
    %end
    %U = Uk;
    licznik = licznik+1;
    error(licznik) = max(max(abs(Uk-U)));

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
%toc