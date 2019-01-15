function [ r ] = createR( n )
a = sin(pi/(2*n))^2;
%disp(a)
k = 1/2 * log(a)/log(sqrt(2)-1);
%disp(k)
k = floor(k)+1;
%disp(k)
r = zeros(1,k);
for i=1:k
    %c = a^((-i)/(k));
    r(i) = a^((1-2*i)/(2*k));
end
end

