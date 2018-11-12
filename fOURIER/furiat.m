clear, clc;
a = 0;
b = 1;

a0 = 2/3;
T = (b-a);
%n = 1;
fan = @(n) 1/(pi^2 * n^2);
%An = fan(n);
fbn = @(n) -1/(pi*n);
%Bn = fbn(n);
%f = @(x) An * cos(2*n*pi/T * x) + Bn * sin(2*n*pi/T * x);
x = linspace(a,b,1000);
S = zeros(size(x),1);
for i=1:5000
  n = i;
  An = fan(n);
  Bn = fbn(n);
  f = @(x) An * cos(2*n*pi/T * x) + Bn * sin(2*n*pi/T * x);
  S += f(x);
  
  
end
%n=1;
%Bn = 1;
%f = @(x) An * cos(2*n*pi/T * x) + Bn * sin(2*n*pi/T * x);
%S += f(x);
S += a0/2;
plot(x,S);
