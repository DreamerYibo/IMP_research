%show the trajectory of the eigenvalue with the largest real part when lambda is changing.
clc; clear;
% Generate (A_i, B_i)

C_1 = [eye(2,2), zeros(2,2)];
C_2 =  [zeros(2,2), eye(2,2)];
S_omega = [0 -1; 1,0];
for (i=1)
    while(true)
        while(true)
            A_11_tilde = diag(4*rand(2,1)); %unstable
            temp_S = 10*(rand(2,2)-0.5);
            A_11_tilde = inv(temp_S)*A_11_tilde*temp_S;
            B_1_tilde = 2*(rand(2,2)-0.5);
            if(rank(ctrb(A_11_tilde, B_1_tilde),1e-3)==2) % controllable
                break;
            end
        end
        
        A_22_tilde = diag(-6*rand(2,1));
        temp_S = 10*(rand(2,2)-0.5);
        A_22_tilde = inv(temp_S)*A_22_tilde*temp_S;
        A_12_tilde = 3*rand(2,2);
        
        P = 6*(rand(4,4)-0.5);
        A = [A_11_tilde, A_12_tilde; zeros(2,2), A_22_tilde];
        A = inv(P)*A*P;
        B = inv(P)*[B_1_tilde;zeros(2,2)];
        
        A=round(A,1); B= round(B,1);

        % %Try this
        % A(:,1:2) = 0;
        if(Is_stabilizable(A,B, 1e-3) && Is_detectable(A,C_1, 1e-3) && IMP_solvability(A,B,C_1,S_omega))
            break
        end
    end
end


setlmis([]);
K1_var=lmivar(2,[2,2]);
K2_var=lmivar(2,[2,2]);
% P = lmivar(1,[4,4]);

lmiterm([1 1 1 K1_var],B,C_1,'s');               % LMI #1: B*K1*C_1 
lmiterm([1 1 1 K2_var],B,C_2,'s');               % LMI #1: B*K2*C_2
lmiterm([1 1 1 0],A'+A);                        % LMI #1: A'+A

SDA=getlmis;

[tmin,xfeas] = feasp(SDA)
K1 = dec2mat(SDA,xfeas,K1_var)
K2 = dec2mat(SDA,xfeas,K2_var)

A_c1 = A+B*K1*C_1+B*K2*C_2;
e_1 = eig(A_c1)
eig(A_c1+A_c1')


N = 1e2;
lambda_list = linspace(1,50, N);
max_eig_real_list= zeros(1,N);

for i=1:N
        lambda = lambda_list(i);
        temp = real(eig(A+lambda*B*K1*C_1+B*K2*C_2));
        max_eig_real_list(i) = max(temp);
end
figure(1)
plot(max_eig_real_list,'x')

%%Part 2

A2 = -50+100*rand(10,10);

N = 1e2;
lambda_list = linspace(1,50, N);
max_eig_real_list= zeros(1,N);

scale_ind = [1:3];

for i=1:N
    lambda = lambda_list(i);
    A2_temp = A2;
    A2_temp(:,scale_ind) = A2(:,scale_ind)*lambda;
    temp = real(eig(A2_temp));
    max_eig_real_list(i) = max(temp);
end

figure(2)
plot(max_eig_real_list,'x')