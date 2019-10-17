function fun = function_from_vector(x, y, arg)

% Younes, 02/02/2016


n=length(arg);

fun = zeros(n,1);

for i= 1:n


z=arg(i)-x;

[~,ind]=min(abs(z));

fun(i) = y(ind);

end
