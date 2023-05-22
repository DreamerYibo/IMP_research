function [PI,GAMMA] = mod_IMP_full_solver(A,B,C,P,Q,S)

    %% For NON SQUARE combined matrix
    %See Chapter 1 of "Topics of Control theory"
    %%Judge if there is a unique solution.
    %See Corollary A1.2. This part is yet finished. Try getting the result at first
    %%Get the result
    n = size(A,1) % size of x
    m = size(B,2) % size of u
    r = size(P,2) %size of omega
    p = size(C,1) %size of e

    PI = zeros(n,r);
    GAMMA = zeros(m,r);

    % if (m~=p)
    %     error("m!=p. There will be more than 1 solution!")
    % end

    A1 = [A, B; C, zeros(p,m)]
    temp = zeros(n+p,n+m);
    temp(1:n, 1:n) = eye(n,n);
    A2 = temp;

    temp2 = kron(transpose(eye(r,r)), A1) -  kron(transpose(S), A2)
    R = [-P;-Q];

    vec_merge = linsolve(temp2, reshape(R,[(n+p)*r,1]))
    matrix_merge = reshape(vec_merge, [n+m, r])
    PI = matrix_merge(1:n, :);
    GAMMA = matrix_merge((n+1):end, :);
end
    
    