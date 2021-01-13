function C = sub_intersect(A,B)
    kernel = null([A -B]); % find null space of [A -B]
    [n,l] = size(A); % detemine dimension of A 
    lambda = kernel(1:l,:); % determine lambda
    C = A*lambda; % x = A*lambda
end

