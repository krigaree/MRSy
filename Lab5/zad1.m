clear, clc;
% zad 1
nx = 50;
ny = 50;
xa = 0;
xb = 1;
ya = 0;
yb = 1;
h = (xb-xa)/(nx+1);
k = (yb-ya)/(ny+1);
xM = linspace(xa + h, xb - h, nx);
yM  = linspace(ya+k, yb-k, ny);

f = @(x, y) 0;
g = @(x, y) x.*y;
uxa = @(x) 0; %u1
uxb = @(x) x; %u3
uya = @(y) 0; %u4
uyb = @(y) y; %u2

T = -4* eye(nx) + diag(diag(eye(nx-1)),-1) + diag(diag(eye(nx-1)),1);
I = eye(nx);
B = kron(eye(ny), T);
C = kron(diag(diag(eye(ny-1)),-1), I);
D = kron(diag(diag(eye(ny-1)),1), I);
A = B + C + D;

[X1 Y1] = meshgrid(xM,yM);
F = f(X1, Y1)' .* ones(nx, ny);
F = h^2 .* F;
F(:, 1) = F(:, 1) - uxa(xM)';
F(:, ny) = F(:, ny) - uxb(xM)';
F(1, :) = F(1, :) - uya(yM);
F(nx, :) = F(nx, :) - uyb(yM);
F = reshape(F,  nx*ny, 1);
U = linsolve(A,F);

U = reshape(U,  nx, ny)';
[X,Y] = meshgrid(xa:h:xb, ya:k:yb);
U = [uya(yM)' .* diag(eye(ny)), U, uyb(yM)' .* diag(eye(ny))];
XM = [xa xM xb];
U = [uxa(XM)' .* diag(eye(nx+2)), U', uxb(XM)' .* diag(eye(nx+2))]';

[X2,Y2] = meshgrid(linspace(xa,xb,100), linspace(ya,yb,100));
G = g(X2,Y2);

Error = max(max(abs(g(X,Y)-U)))
subplot(1,2,1)
surf(X,Y,U)
subplot(1,2,2)
surf(X2,Y2,G)