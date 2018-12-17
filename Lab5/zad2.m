clear, clc;
% zad 3
% Przykładowe dane
nx = 10; ny = 10; xa = 0; xb = 2; ya = 0; yb = 1;
f = @(x, y) (x.^2 + y.^2) .* exp(x.*y);
g = @(x, y) exp(x.*y);
uxa = @(x) 1;         %u1
uxb = @(x) exp(x);    %u3
uya = @(y) 1;         %u4
uyb = @(y) exp(2.*y); %u2

% Obliczenia
h = (xb-xa)/(nx+1);
k = (yb-ya)/(ny+1);
xM = linspace(xa + h, xb - h, nx);
yM  = linspace(ya+k, yb-k, ny);
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
% Rozwiązanie
U = linsolve(A,F);
U = reshape(U,  nx, ny)';
[X,Y] = meshgrid(xa:h:xb, ya:k:yb);
U = [uya(yM)' .* diag(eye(ny)), U, uyb(yM)' .* diag(eye(ny))];
XM = [xa xM xb];
G = g(X,Y);
U = [uxa(XM)' .* diag(eye(nx+2)), U', uxb(XM)' .* diag(eye(nx+2))]';
% Błąd
Error = max(max(abs(G-U)));
% Wykres
subplot(1,2,1)
surf(X,Y,U)
subplot(1,2,2)
surf(X,Y,G)
