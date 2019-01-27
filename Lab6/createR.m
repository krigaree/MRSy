function [ r ] = createR( n )
a = sin(pi/(2*n))^2;
k = 1/2 * log(a)/log(sqrt(2)-1);
k = floor(k)+1;
r = zeros(1,k);
for i=1:k
    r(i) = a^((1-i)/(2*k));
end
end

