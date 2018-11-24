clear, clc;
#dane wej�ciowe
a=0;
b=1;
n=input('Podaj liczbe wezlow: ');
h = (b-a)/(n-1);
Ua = 0;
Ub = 1;
f = @(x) 12*x;
g = @(x) 2*x.^3-x;
#obliczenia
v1 = -2*diag(eye(n));
v2 = diag(eye(n-1));
A = diag(v1) + diag(v2,1) + diag(v2,-1);
x= linspace(a+h,b-h,n);
x2 = linspace(a,b,n+2);
F = f(x)*h^2;
F(1) = F(1) - Ua;
F(n) = F(n) - Ub;
U = linsolve(A,F');
U = [Ua U' Ub];
#wykres
hold on
plot(x2,U', 'r')
plot(x2,g(x2),'b')
hold off