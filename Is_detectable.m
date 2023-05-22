function true_or_false = Is_detectable(A,C, tol)
    Lambda_list = eig(A);
    Lambda_list = Lambda_list(real(Lambda_list)>0 | real(Lambda_list)==0);
    N = length(Lambda_list);
    if(N == 0)
        true_or_false = 1;
    end
    n = size(A,1); % size of x
    % m = size(B,2); % size of u
    for i=1:N 
        M = [A-Lambda_list(i)*eye(n); C];
        if (rank(M, tol) ~= n) %tolerance of rank(). Eg, 1e-4
            true_or_false = 0;
        return
    end
    true_or_false = 1;

end