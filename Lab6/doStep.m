function [Uv] = doStep (iter, r, f, x, y, U, step)
if step %n-step -> n+1/2
  m = length(x);
  fi = U(2:m-1, iter) + r .* (U(2:m-1,iter+1) + U(2:m-1,iter-1) - 2*U(2:m-1,iter));
  F = fi .* diag(eye(m-2));
  F(1) = F(1) + r .* U(1, iter);
  F(m-2) = F(m-2) + r .* U(m, iter);
  U(m-2, iter);
else %n+1/2-step -> n
  m = length(y);
  fi = U(iter, 2:m-1) + r .* (U(iter+1, 2:m-1) + U(iter-1, 2:m-1) - 2*U(iter, 2:m-1));
  F = fi' .* diag(eye(m-2));
  F(1) = F(1) + r .* U(iter, 1);
  F(m-2) = F(m-2) + r .* U(iter, m);
  U(iter, m-2);
end
m = m-2; %delete x0 xm+1
A = (1+2*r).*eye(m) + diag(-r*diag(eye(m-1)),-1) + diag(-r*diag(eye(m-1)),1);
%disp(size(A))
A;
F;
Uv = linsolve(A,F);
  

end
