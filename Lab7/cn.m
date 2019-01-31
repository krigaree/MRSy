clc,clear all
tic
%rozwi¹zanie analityczne
G = @(x,t) sin(pi.*x./2).*exp(-(pi.^2).*t./4);

%przedzia³ omega
xa=0;
xb=2;
yc=0;
yd=1;

%warunki brzegowe
ub1 = @(t) 0;
ub2 = @(t) 0;
up = @(x) sin(pi*x/2);

%siatka
n=50;

h=(xb-xa)/(n+1);
k=(yd-yc)/(n+1);
x=xa:h:xb;    
t=yc:k:yd;

D=1;
alfa = (D*k)/(2*(h^2));

%tworzenie macierzy A
A1 = -alfa*ones(n+1,1);
A2 = (2*alfa+1)*ones(n+2,1);
A = diag(A1,-1) + diag(A2) + diag(A1,1);
A(1,1)=1; A(1,2)=0;
A(end,end)=1; A(end,end-1)=0;

%utworzenie pustej macierzy
U=zeros(n+2,n+2); 

%dodanie warunków pocz¹tkowego/brzegowego
U(1,:) = up(x);
U(:,1) = ub1(t);
U(:,n+2) = ub2(t);

psi = zeros(1,n+2);

licznik=0;

for j=2:n+2
    for i=2:n+1
        psi(i) = alfa*(U(j-1,i+1)-2*U(j-1,i)+U(j-1,i-1))+U(j-1,i);
    end
    
    F = linsolve(A,psi');
    U(j,1:end)=F;
    licznik=licznik+1;
end

[X,T] = meshgrid(x,t);
subplot(1,2,1)
surf(X,T,U)
title('Metoda Numeryczna')
subplot(1,2,2)
surf(X,T,(G(X,T)))
title('Metoda Analityczna')
Error=max(max(abs(U-G(X,T))))
licznik
toc