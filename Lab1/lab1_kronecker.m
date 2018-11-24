clc,clear;

A=[2, -1;-1,2];
B=[0,-1;2,3];
C=[1,1;1,2];

n=6;

%E1
E1 = kron(eye(n),A);

%E2
v = [1;2*diag(eye(n-2));1];

E2 = kron(diag(v),A);

%E3
v = [1;-1*diag(eye(n-2));1];

E3 = kron(diag(v),A) +\
  kron(diag(diag(eye(n-1)),1),B) +\
  kron(diag(diag(eye(n-1)),-1),B);
  
%E4  
v = [1;2*diag(eye(n-2));1];

E4 = kron(diag(v),A) +\
  kron(diag(diag(eye(n-1)),1),B) +\
  kron(diag(diag(eye(n-1)),-1),B');

%E5
E5=kron(eye(n),A) +\
  kron(diag(diag(eye(n-2)),2),B') +\
  kron(diag(diag(eye(n-2)),-2),B);

%E6  
v = [1;zeros(n-2,1);1];
v1 = [0;diag(eye(n-2));0];

E6 = kron(diag(v),C)+ kron(diag(v1),A) +\
  kron(diag(diag(eye(n-2)),2),3*B') +\
  kron(diag(diag(eye(n-2)),-2),2*B) +\
  flip(kron(diag(v), flip(eye(2))));
  
%E7
v=[1;zeros(n-2,1);1];
v1=[0;1;zeros(n-4,1);1;0];
v2=[0;0;2*diag(eye(n-4));0;0];
 
E7=kron(diag(v),A) +\
  kron(diag(v1),B) +\
  kron(diag(v2),B) +\
  kron(diag(diag(eye(n-2)),2),C) +\
  kron(diag(diag(eye(n-2)),-2),C');
  
%E8
v=[1;zeros(n-2,1);1];
v1=[0;2;zeros(n-4,1);2;0];
v2=[0;0;diag(eye(n-4));0;0];

E8=kron(diag(v),A) +\
  kron(diag(v1),B) + kron(diag(v2),C) +\
  kron(diag(diag(eye(n-2)),2),eye(2)) +\
  kron(diag(diag(eye(n-2)),-2),eye(2))+\
  flip(kron(diag(v),flip(A)));
  
E1
E2
E3
E4
E5
E6
E7
E8