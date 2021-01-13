function C = sub_sum(A,B)
    M = [A B]; % V + W = Im([V W])
    C = orth(M); % return basis for M
end

