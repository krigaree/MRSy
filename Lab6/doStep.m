function [Uv] = doStep (iter, r, x, y, U, step)
if step %n-step -> n+1/2  
  m = length(x);
  fi = U(iter, 2:m-1) + r .* (U(iter+1, 2:m-1) + U(iter-1, 2:m-1) - 2*U(iter, 2:m-1));
  F = fi' .* diag(eye(m-2));
  F(1) = F(1) + r .* U(iter, 1);
  F(length(F)) = F(length(F)) + r .* U(iter, m);
else %n+1/2-step -> n   
  m = length(y);
  fi = U(2:m-1, iter) + r .* (U(2:m-1,iter+1) + U(2:m-1,iter-1) - 2*U(2:m-1,iter));
  F = fi .* diag(eye(m-2));
  F(1) = F(1) + r .* U(1, iter);
  F(length(F)) = F(length(F)) + r .* U(m, iter);
end
m = m-2; %delete x0 xm+1
A = (1+2*r).*eye(m) + diag(-r*diag(eye(m-1)),-1) + diag(-r*diag(eye(m-1)),1);
Uv = linsolve(A,F);
end
