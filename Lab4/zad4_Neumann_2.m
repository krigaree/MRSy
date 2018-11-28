clear, clc;
f = @(x) exp(x) .* (-x .^2 + x + 2);
g = @(x) x .* exp(x);
a = 0;
b = 1;
ax = @(x) 1;
bx = @(x) -x;
cx = @(x) 1;
ua = 0;
ubf = 2 * exp(1);
c = b-a;
n = 25;
h = c/(n+1);
x = linspace((a+h),(b-h),n);

a1 = [(-2 * ax(x) + cx(x) * h^2).* diag(eye(n)); -1];
A1 = diag(a1);
a2 = (ax(x) +1/2 * bx(x) * h)' .* diag(eye(n));
A2 = diag(a2, 1);
a3 = [(ax(x(2:n)) -1/2 * bx(x(2:n)) * h)' .* diag(eye(n-1)); 1];
A3 = diag(a3, -1);

A = A1  + A2 + A3;

F = h^2 * f(x);
F(1) = F(1) - ua;
F = [F 1/2 * h^2*f(b) - h * ubf];
U = linsolve(A,F');
U = [ua U'];
X = [a x b];
plot(X, U, X, g(X), 'ro');
E = max(abs(g(X) - U));