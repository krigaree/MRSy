clear, clc;
#dane wejœciowe
a=0;
b=2*pi;
n=input('Podaj liczbe wezlow: ');
h = (b-a)/(n-1);
Ua = 0;
Ub = 0;
f = @(x) -sin(x) - 4*sin(2*x);
g = @(x) sin(x) + sin(2*x);
#obliczenia
v1 = -2*diag(eye(n));
v2 = diag(eye(n-1));
A = diag(v1) + diag(v2,1) + diag(v2,-1);
x= linspace(a+h,b-h,n);
x2 = linspace(a,b,n+2);
x3 = linspace(a,b,100);
F = f(x)*h^2;
F(1) = F(1) - Ua;
F(n) = F(n) - Ub;
U = linsolve(A,F');
U = [Ua U' Ub];
#wykres
plot(x3, g(x3), x2, g(x2), 'ro')
legend('Metoda Analityczna', 'Metoda Numeryczna');
#error
E = max(abs(g(x2) - U))
