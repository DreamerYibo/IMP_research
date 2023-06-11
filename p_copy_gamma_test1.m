% Add one column in B. Adjust L such that Gamma_1 is 0.
clc;clear all;

map_size = [100,100]; % cm
N=1;n=4; %agent 1 original L_bar, agent i L_bar + delta_L_i
w=1;
l = 7;
S = [0, 1; 0, 0]; % exp(S_t) is counter clockwise ()
Q = [1 -1; 1 1]; % y1 ---> Qv. y ---> [t; t]
v0 = [0 4];


A_bar = [0 1 0 0;...
         0 0 0 0;...
         0 0 0 1;...
         0 0 0 0];

% A = A_bar;

% A = A_bar + rand(n,n)-0.5;
% A = A_bar -[ 0   -0.1546    0.1973    0.2723;
%             -0   -0.1459   -0.1236    0.2417;
%             -0    0.0647    0.2809   -0.1941;
%             -0  -0.0933    0.1338   -0.1169];

            % A = A_bar +[  0     1     0.1973      -0.029077;
            %             0     0     -0.1236    -0.29226;  
            %             0     0     0.2809     -0.1941;
            %             0     0     0.1338     -0.1169;];

% A = A_bar + [0.5-rand(4,2),zeros(4,1),0.5-rand(4,1)];
A = A_bar +  [0   -0.2749         0    0.1002;
            0  -0.1         0    0.1;
            -0   -0.1        0   -0.1;
            0    0         0    0;];


% A = -A_bar;

% B = [0 0;
%      1 0;
%      0 0; 
%      0 1] + (0.5-rand(4,2));
% B_bar = [0 1;
%         1 0;
%         1 0; 
%         0 1];

B_bar = [0 0 0;
        1 0 1;
        0 0 1; 
        0 1 1];
B = B_bar;
% B = B_bar + [ 0.0073   -0.2079;
%      -0.1501   -0.1461;
%       0.1791   -0.2297;
%       0.0391   -0.0779];

% B = B_bar + [ 0.0073   -0.2079;
%      -0.1501   -0.1461;
%       0.1791   -0.2297;
%       0.0391   -0.0779];

B = B_bar+(0.5-rand(4,3))*0.1;

C_bar= [1 0 0 0;
        0 0 1 0];

% C = C_bar;
C = C_bar +0.5*(rand(2,4)-0.5);



A_c = kron(eye(N), A);
B_c = kron(eye(N), B);
C_c = kron(eye(N), C);

% p-copy
% n = size(A,1); % size of x
m = size(B,2); % size of u
r = size(S,2); %size of omega
p = size(C,1); %size of e

G_1 = kron(eye(p), [0 1; 0 0]);
G_2 = kron(eye(p), [0; 1]);

[~, K, ~] = icare([A_bar, zeros(n,r*p); G_2*C_bar, G_1], [B_bar; zeros(r*p,m)] ,3*eye(n+r*p),[],[],[],[]); %random stabilizing sol
K = -K;



L = zeros(m,r); % set as 0 at first
%%some other tests: start

A1 = [A, B; C, zeros(p,m)];
    temp = zeros(n+p,n+m);
    temp(1:n, 1:n) = eye(n,n);
    A2 = temp;

    temp2 = kron(transpose(eye(r,r)), A1) -  kron(transpose(S), A2);
    temp2_b = reshape([zeros(n,r); Q],[(n+p)*r,1]);
    null_temp2 = null(temp2);
    null_pi_gamma = reshape(null_temp2, [n+m,r,size(null_temp2, 2)])
    pi_gamma_vec=pinv([temp2; [zeros(m,n),eye(m),zeros(m,n+m)]])*[temp2_b; zeros(m,1)]; % force the first col of Gamma to be zero.

    if (norm([temp2; [zeros(m,n),eye(m),zeros(m,n+m)]]*pi_gamma_vec - [temp2_b; zeros(m,1)])>0.01)
        [temp2; [zeros(m,n),eye(m),zeros(m,n+m)]]*pi_gamma_vec - [temp2_b; zeros(m,1)]
        error('cannot have a correct solution for Gamma') 
    end

   

    pi_gamma = reshape(pi_gamma_vec, [n+m,r]);
    Pi = pi_gamma(1:n,:)
    Gamma = pi_gamma(n+1:end, :)

%%some other tests: end

% mapping from l to gamma_1 based on actual model (see "B*L as the input)
A_f_actual =[A, zeros(n,r*p); G_2*C, G_1] + [B; zeros(r*p,m)]*K;
actual_sylv_sq_mat = -kron(S', eye(2*n,2*n))+kron(eye(2,2)', A_f_actual);
actual_px_pl = - inv(actual_sylv_sq_mat)*kron(eye(r), [B; zeros(r*p,m)]);
actual_px_pl_offset =  - inv(actual_sylv_sq_mat)*reshape([zeros(n,2);-G_2*Q], [r*n+r*p*r, 1]);

actual_pgamma_pl = kron(eye(r),K)*actual_px_pl + eye(m*r,m*r);
actual_pgamma_pl_offset = kron(eye(r),K)*actual_px_pl_offset;

actual_pz_pl = [eye(m), zeros(m,m)]*actual_pgamma_pl;
actual_pz_pl_offset = [eye(m), zeros(m,m)]*actual_pgamma_pl_offset;

actual_gamma = actual_pz_pl*L(:)+actual_pz_pl_offset

% mapping from l to gamma_1 based on nominal model (see "B*L as the input)
A_f_nominal =[A_bar, zeros(n,r*p); G_2*C_bar, G_1] + [B_bar; zeros(r*p,m)]*K;
nominal_sylv_sq_mat = -kron(S', eye(2*n,2*n))+kron(eye(2,2)', A_f_nominal);
nominal_px_pl = - inv(actual_sylv_sq_mat)*kron(eye(r), [B_bar; zeros(r*p,m)]); %from l to pi_sigma
nominal_px_pl_offset =  - inv(actual_sylv_sq_mat)*reshape([zeros(n,2);-G_2*Q], [r*n+r*p*r, 1]);

nominal_pgamma_pl = kron(eye(r),K)*nominal_px_pl + eye(m*r,m*r);
nominal_pgamma_pl_offset = kron(eye(r),K)*nominal_px_pl_offset;

nominal_pz_pl = [eye(m), zeros(m,m)]*nominal_pgamma_pl;
nominal_pz_pl_offset = [eye(m), zeros(m,m)]*nominal_pgamma_pl_offset;

pz_pl = nominal_pz_pl;
pz_pl_offset = nominal_pz_pl_offset;




dt = 0.01;
t_max = 300;
t=0

X = zeros(N,t_max/dt); % x position
Y = zeros(N,t_max/dt);
% leader_offset = [30; 30];
Omega_c = zeros(2*N,t_max/dt);
X_c =  zeros(n*N,t_max/dt);
Xi_c = zeros(r*p*N,t_max/dt);

Gramm_c = [];% obtained from [t, t+1]
epsilon_list = [];
epsilon_time = [];
R_relative = zeros(N,t_max/dt);

Omega_c(:,1) = v0;

X(:,1) =0.3*map_size(1)*(rand(N,1)-0.5);
Y(:,1) =0.3*map_size(2)*(rand(N,1)-0.5);

% disp("Initial  positions")
% [X(:,1)'; Y(:,1)']


X_c(1:n:end, 1) = X(:,1); %All other states are set to zero. Xi_c also set to all zero
X_c(2:n:end, 1) = 0;
X_c(3:n:end, 1) = Y(:,1); % for double integrator like A B C
X_c(4:n:end, 1) = 0;

disp("Initial conditions")
reshape(X_c(:,1),[n,N])


A_closed = [[A, zeros(n,r*p); G_2*C, G_1] + [B; zeros(r*p,m)]*K, [B*L; G_2*(-Q)];
            zeros(r, n+r*p), kron(eye(N), S);];

state_cl = [X_c(:,1); Xi_c(:,1); Omega_c(:,1)]; 

% Prof's idea, optimization like
h_gram = 1; %grammian calculation period
counter = 0; %count the times of grammian calculations
update_L_index = 0;

%For updating the map from l to z
l_list = [];
z_list = [];

d_L = 0;

for k= 2:(t_max/dt)
    t = t+dt;
    state_cl = expm(A_closed*dt)*state_cl;

    X_c(:,k) = state_cl(1:n);
    Xi_c(:,k) = state_cl(n+1:n+r*p);
    Omega_c(:,k) = state_cl(n+r*p+1:end);
        
    for j =1:N
        X(j,k) = X_c(n*(j-1)+1, k); % Assume C = (I_2 0)
        Y(j,k) = X_c(n*(j-1)+3, k);
    end
    if (mod(t,h_gram) < dt) %calculate grammian if (mod(t,1) == 0) does not work, very strange.
        counter = counter+1; %how many times gramm is calculated
        k2 = k;
        k1 = k-h_gram/dt;
        
        temp_gramm = zeros(m*N,2);
        
        
        for j=(1:N)
            itgl_1 = zeros(r,r); %integral
            itgl_2 = zeros(m,r);
            for l=(k1:k2-1)
                u = K*[X_c(n*(j-1)+1: n*j,l);Xi_c(n*(j-1)+1: n*j,l)]+L*Omega_c(2*j-1:2*j,l);
                itgl_2 = itgl_2 + u*Omega_c(2*j-1:2*j,l)'*dt;
                itgl_1 = itgl_1 + Omega_c(2*j-1:2*j,l)* Omega_c(2*j-1:2*j,l)'*dt;
            end
            temp_gramm(m*(j-1)+1:m*j, :) = itgl_2*inv(itgl_1);
        end
        
        
        Gramm_c = [Gramm_c, temp_gramm];
        
        if (counter > 1 && counter - update_L_index > 1)
            
            
            if (all( abs(Gramm_c(:,end-1:end) - Gramm_c(:,end-3:end-2)) < 0.01) )
                % % norm(Gramm_c(:,end-1:end))*0.001
                % %debug: get the pi from the actual model
                % pi_sigma = -inv(actual_sylv_sq_mat)*kron(eye(2,2), blkdiag(B, eye(n,n)))*reshape([L_1_bar; L_2_bar], [2*(n+2),1]);
                % temp_gramm = reshape(pi_sigma, [2*n, 2]);
                % temp_gramm = temp_gramm(1:n,:);
                
                % temp_gramm - Gramm_c(:,end-1:end);
                % Gramm_c(:,end-1:end) = temp_gramm;
                % %debug end
                
                l_temp = reshape(L, [m*r,1]);
                
                if (update_L_index == 0 || true )
                    if (norm(d_L)/norm(L) > 0.01)
                        disp("pz_pl updated, time: "+ num2str(t))
                        l_list = [l_list, l_temp];
                        z_list = [z_list, reshape(Gramm_c(:,end-1), [m,1])]; %only get first col of Gamma
                        
                        temp_mat = kron([l_list; ones(1, size(l_list, 2))]', eye(m,m));
                        nominal_R_bar = [nominal_pz_pl, nominal_pz_pl_offset]; %estimate bothe pz_pl and its offset
                        delta_R_bar_update = pinv(temp_mat)*(z_list(:) - temp_mat*nominal_R_bar(:)) % get the solution closest to the nominal one when there are many possible solutions. When the size gets larger, the one get the linear least square sol.
                    
                        pz_pl = reshape(nominal_pz_pl(:)+delta_R_bar_update(1:end - m), [m,m*r]);
                        pz_pl_offset = reshape(nominal_pz_pl_offset(:)+delta_R_bar_update(end-m+1:end), [m,1]);
                    end
                    z_1 = reshape(Gramm_c(:,end-1), [m,1]); % z = vec(C*Pi)
                    l_1 = reshape(L, [m*r,1]);
                    
                    % if (rank([l_list; ones(1, size(l_list, 2))]') == m*r+1)
                    if (rank([l_list; ones(1, size(l_list, 2))]') == m*r+1 && norm(d_L)/norm(L) < 0.01)
                        disp("I am here")
                        % d_l_2 = pinv(pz_pl)*(0 - pz_pl_offset-pz_pl*l_1); %closest one to the previous l
                        d_l = pinv(pz_pl)*(-z_1-pz_pl*l_1);
                        l_2 = l_1+d_l;
                    else
                        step_length = 1;
                        
                        pe_pz = z_1';
                        pe_pl = pe_pz*pz_pl;
                        delta_l = -pe_pl';
                        
                        l_2 = l_1 + step_length*delta_l;
                        % Pi_2 = reshape(reshape(Pi_1,[2*n,1]) + step_length*delta_pi, [n,2] );
                        z_2 = z_1 + pz_pl*(l_2-l_1); %estimated change in z, map: z->l is linear and it is pz_pl (based on nominal model).
                        
                        e_1 = 1/2*(z_1'*z_1);
                        e_2 = 1/2*(z_2'*z_2);
                        
                        epsilon_list = [epsilon_list, e_1];
                        epsilon_time = [epsilon_time, t];
                        
                        % e_2-e_1
                        while (e_2-e_1 > 0.3 * step_length* pe_pl * delta_l)
                            step_length = step_length/2;
                            
                            l_2 = l_1 + step_length*delta_l;
                            z_2 = z_1 + pz_pl*(l_2-l_1);
                            
                            e_2 = 1/2*(z_2'*z_2);
                            % if (step_length < 0.01)
                            %     break;
                            % end
                        end
                    end
                    
                    %update the L
                    
                    d_L = reshape(l_2-l_1, [m,2])
                    L = L + d_L;
                    % L_2_bar = L_2_bar + d_L(3:end,:);
                    
                    disp("L1 L2 updated, time: "+ num2str(t))
                    % L_bar = L_bar + d_L;
                    
                    
                    %DEBUG
                    % Pi_2
                    % temp_sol = linsolve(kron(S', eye(n,n))-kron(eye(2,2)', A+B*K), reshape(B*L_bar, [n*2,1]));
                    
                    % Pi_2_actual = reshape(temp_sol, [n,2])
                    % e_CPi_2_actual = norm(C*Pi_2_actual-eye(2), 'fro')^2;
                    
                    %update closed-loop sys
                    % L_1c = kron(eye(N),L);
                    % L_2c = kron(eye(N),L_2_bar);
                    
                    A_closed = [[A, zeros(n,r*p); G_2*C, G_1] + [B; zeros(r*p,m)]*K,        [B*L;      G_2*(-Q)];
                    zeros(r, n+r*p), kron(eye(N), S);];
                    
                    
                    
                    update_L_index = counter;
                    
                    % W_op = diag([t,t,1,t]);
                    
                    % % debug
                    % z_temp = actual_pz_pl*l_temp;
                    % z_temp - reshape(C*Gramm_c(:,end-1:end), [4,1])
                    
                    % disp("error of L regression")
                    % norm(z_list(:) - temp_mat*pz_pl(:))
                    % disp("actual error")
                    % norm(z_list(:) - temp_mat*actual_pz_pl(:))
                    
                end
            end
        end
    end
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end

% rgb_list = rand(N,3);
rgb_list = [255 128 0; 0 153 0; 0 204 204; 204 204 0] * 1/255;

E_norm = zeros(N,t_max/dt);
U_norm = zeros(N,t_max/dt);
for k = 1:t_max/dt
    % u = K*[X_c(:,k);Xi_c(:,k)]+L*Omega_c(:,k);
    for i = 1:N 
        E_norm(i,k) = norm(C*X_c((i-1)*n + 1:  i*n, k) - Q*Omega_c(2*i-1: 2*i ,k));
        % U_norm(i,k) = norm(u((i-1)*m+1: i*m));
    end
end
figure(5)
clf(figure(5))
% E = C_c*X_c+Q_c*Omega_; % without using a leader. So the center of circle can be somewhere related to the initial positions of agents
for (i=1:N)
    temp_h(i) = plot(dt:dt:t_max , E_norm(i,:), 'Color', rgb_list(i,:), 'LineWidth', 1);
    hold on
    str_array(i) = "Agent"+num2str(i);
end
legend(temp_h, num2cell(str_array));
xlabel('time $t \,\, (s)$', 'interpreter', 'latex')
ylabel("$\parallel y_i-v_i-\bar{v} \,\,\parallel$", 'interpreter', 'latex')

% figure(4)
% clf(figure(4))
% for (i=1:N)
%     temp_h(i) = plot(dt:dt:t_max , U_norm(i,:), 'Color', rgb_list(i,:), 'LineWidth', 1);
%     hold on
%     str_array(i) = "Agent"+num2str(i);
% end
% legend(temp_h, num2cell(str_array));
% xlabel('time $t \,\, (s)$', 'interpreter', 'latex')
% ylabel("$\parallel u_i \,\,\parallel$", 'interpreter', 'latex')

figure(6)
clf(figure(6))
% E = C_c*X_c+Q_c*Omega_; % without using a leader. So the center of circle can be somewhere related to the initial positions of agents
plot(epsilon_time , epsilon_list, '-x','Color', rgb_list(1,:), 'LineWidth', 1);
