clear, clc;
f = @(x) -4*x.^2;
g = @(x) x.^2 - x;
a = 0;
b = 1;
ax = @(x) -x.^2;
bx = @(x) -2*x;
cx = @(x) 2;
ua = 0;
ub = 0;
c = b-a;
n = 10;
h = c/(n+1);
x = linspace((a+h),(b-h),n);

a1 = (-2 .* ax(x) + cx(x) .*h^2)' .* diag(eye(n));
A1 = diag(a1);
a2 = (ax(x(2:n)) -1/2 .* bx(x(2:n)) .*h)' .* diag(eye(n-1));
A2 = diag(a2, -1);
a3 = (ax(x(1:n-1)) +1/2 .* bx(x(1:n-1)) .*h)' .* diag(eye(n-1));
A3 = diag(a3, 1);

A = A1  + A2 + A3;

F = (h^2 * f(x));
F(1) = F(1) - ua * (ax(1) - 1/2 * bx(1)*h);
F(n) = F(n) - ub * (ax(n) + 1/2 * bx(n)*h);
U = linsolve(A,F');
U = [ua U' ub];
X = [a x b];
plot(X, U, X, g(X), 'ro');
E = max(abs(g(X) - U));