clc; clear;
% Generate (A_i, B_i)
N=4;
A_list = zeros(4,4,N);
B_list = zeros(4,2,N);
K_list = zeros(2,4,N);
C_1 = [eye(2,2), zeros(2,2)];
S_omega = [0 -1; 1,0];
for (i=1:4)
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
            A_list(:,:,i) = A;
            B_list(:,:,i) = B;
            [~,K,~]=icare(A,B,100*[eye(2),zeros(2,2);zeros(2,4)]+eye(4,4),eye(2,2),zeros(4,2),eye(4,4),zeros(4,4));
            % [~,K,~]=icare(A,B,0.1*eye(4,4),eye(2,2),zeros(4,2),eye(4,4),zeros(4,4));
            K_list(:,:,i) = -K;
            break
        end
    end
end
save("ABK_list")