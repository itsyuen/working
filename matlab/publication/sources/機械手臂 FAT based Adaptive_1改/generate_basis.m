function basis = generate_basis(t, l, T)

if mod(l,2) == 0
    Message = MException('generate_basis:excption','argument l should be odd value, l = %d', l);
    throw(Message);
end

basis = zeros(l, length(t));
n = (1:l);
remainder = mod(n,2);
odd = find(remainder);
even = find(~remainder);
k = (1:length(even)).';

basis(1,:) = ones(1, length(t));
basis(odd(2:end), :) = sin((2*pi/T)*k.*t);
basis(even, :) = cos((2*pi/T)*k.*t);
end
