function true_or_false = IMP_solvability(A,B,C,S)
    %See Chapter 1 of "Topics of Control theory"
    %%Judge if there is a unique solution.
    %See Corollary A1.2. 

    Lambda_list = eig(S);
    N = length(Lambda_list);

    n = size(A,1); % size of x
    m = size(B,2); % size of u
    r = size(S,1); %size of omega
    p = size(C,1); %size of e
    
    if (m~=p)
        error("m!=p. There will be more than 1 solution!")
    end

    for i=1:N 
        M = [A-Lambda_list(i), B; C, zeros(p, m)];
        if (rank(M, 1e-4) ~= (n+p))
            true_or_false = 0;
            return
    end
    true_or_false = 1;

end
    
    