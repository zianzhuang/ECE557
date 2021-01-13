%% Part 1: Basis Vector for Subspaces
clear;
clc;
V = [1 1 3; -1 1 1; 0 0 0; 1 0 1];
W = [1 1; 0 0; 2 -2; 1 0];
V_basis = orth(V)
W_basis = orth(W)

%% Part 2: Sum and Intersection of V and W
VW_sum = sub_sum(V_basis, W_basis)
VW_intersect = sub_intersect(V_basis, W_basis)
rank_VW = rank([V_basis W_basis]) % V and W are not linearly independent

%% Part 3: Coordinates of an R2 vector space
R2_basis = [1 2; 1 1]; 
x = [2;1];
coordinates = R2_basis\x

%% Part 4: matrix representation of linear transformation
A = [1 2 0 -1; 0 1 -1 -2; 1 0 3 4; 0 -1 2 3; 0 0 2 2];
P = [1 1 1 1; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Q = [1 1 0 0 1; 1 -1 0 0 0; 0 0 1 1 0; 0 0 1 -1 0; 0 0 0 0 1];
A_hat = Q\(A*P)

%% Part 5: injective or surjective?
rank_A = rank(A) % rank of A
[m,n] = size(A); % n
nullity_A = n - rank_A % nullity of A

% surjectivity injectivity test
surj = 0;
inj = 0;
if rank_A == m % criteria for surjectivity
    surj = 1; 
end
if nullity_A == 0% criteria for injectivity
    inj = 1;
end

% print test result
if surj == 1 && inj == 1
    disp('A is both surjective and injective')
elseif surj == 1 && inj == 0
    disp('A is surjective, not injective')
elseif surj == 0 && inj == 1
    disp('A is injective, not surjective')
else 
    disp('A is not injective, nor surjective')
end

%% Part 6: number of solutions?
b = [1;0;0;0;0];%
if rank([A b]) > rank(A)
    fprintf('for b = [%d;%d;%d;%d;%d], Ax=b has no solution',b);
elseif nullity_A == 0
    solution = A\b;
    fprintf('for b = [%d;%d;%d;%d;%d], Ax=b has one solution: [%.2f;%.2f;%.2f;%.2f]', b, solution);
else 
    null_basis = null(A);
    solution1 = A\b + null_basis;
    solution2 = A\b - null_basis;
    fprintf('for b = [%d;%d;%d;%d;%d], Ax=b has infinite solutions',b);
    fprintf('\n')
    fprintf('exampls: [%.2f;%.2f;%.2f;%.2f], [%.2f;%.2f;%.2f;%.2f]',solution1, solution2);
end

b = [1;1;-2;-2;-2];
if rank([A b]) > rank(A)
    fprintf('for b = [%d;%d;%d;%d;%d], Ax=b has no solution',b);
elseif nullity_A == 0
    solution = A\b;
    fprintf('for b = [%d;%d;%d;%d;%d], Ax=b has one solution: [%.2f;%.2f;%.2f;%.2f]', b, solution);
else 
    null_basis = null(A);
    solution1 = A\b + null_basis;
    solution2 = A\b - null_basis;
    fprintf('for b = [%d;%d;%d;%d;%d], Ax=b has infinite solutions',b);
    fprintf('\n')
    fprintf('exampls: [%.2f;%.2f;%.2f;%.2f], [%.2f;%.2f;%.2f;%.2f]',solution1, solution2);
end

%% Part 7: A-invariant and Representation Theorem
A = [1 2 2; 1 -3 -5; -1 2 4];
v = [0 1; 1 -2; -1 1];
rank_v = rank(v)
rank_Av_and_v = rank([A*v v])

% Find P for representation theorem
w = null(transpose(v));
P = [v w]
pAP = inv(P)*A*P

%% Part 8: Kalman Decomposition
A = [5 -11 5; 0 -6 0; -5 5 -1];
B = [1 -2; 0 0; 1 2];
Q = ctrb(A,B);
rank_Q = rank(Q)
Q_basis = orth(Q);
w = null(transpose(Q_basis));
P = [Q_basis w]
A_hat = inv(P)*A*P
B_hat = inv(P)*B

% Note: orders are changed for proper "fprintf" representation
ct_1 = [A_hat(1,1) A_hat(2,1);A_hat(1,2) A_hat(2,2)];
ct_2 = [A_hat(1,3);A_hat(2,3)];
ct_3 = [B_hat(1,1) B_hat(2,1); B_hat(1,2) B_hat(2,2)];
unct = A_hat(3,3);

fprintf('The controllable subsystem is z1_dot = [%.2f %.2f;%.2f %.2f]z1+[%.2f;%.2f]z2+[%.2f %.2f; %.2f %.2f]u', ct_1, ct_2, ct_3)
fprintf('\n')
fprintf('The uncontrollable subsystem is z2_dot = %.2f z2',unct)