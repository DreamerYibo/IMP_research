% Unknown A B C dynamics. Output feedback feedback case. Update L1 only. F and S have no common eigenvalues. 
%S is jordan block or unstable

% Change A such that the center can be shifted
clc;clear all;
% load('uncertain_output_test2_bad_result_1', 'A', 'B', 'C')
rot_2d = @(t) [cos(t), -sin(t); sin(t), cos(t)];

map_size = [100,100]; % cm
N=1;n=4; %agent 1 original L_bar, agent i L_bar + delta_L_i
w=1;
l = 7;
S = [0, 1; 0, 0]; %Jordan block
Q = [1 -1; 1 -1]; % y ---> Qv. y ---> [t; t]
v0 = [0 1];

A_bar = [0 1 0 1;...
     0 0 1 -1;...
     0 0 -1 1;...
     0 0 0 0];

% A = A_bar;

% A = A_bar + rand(n,n)-0.5;
% A = A_bar + 0.5*(rand(4,4)-0.5);
% A = A_bar +0.5*[ 0.2786   -0.1546    0.1973    0.2723;
%             -0.2095   -0.1459   -0.1236    0.2417;
%             -0.2604    0.0647    0.2809   -0.1941;
%             -0.1072   -0.0933    0.1338   -0.1169];

A = A_bar +0.5*[ 0   -0.1546    0.1973    0.2723;
            -0   -0.1459   -0.1236    0.2417;
            -0    0.0647    0.2809   -0.1941;
            -0   -0.0933    0.1338   -0.1169];

% A = -0.0247    0.8262   -0.2109    0.7523
%    -0.2081    0.1629    0.9713   -0.8625
%    -0.1355    0.0192   -1.1967    1.1587
%     0.2067    0.2481    0.2309    0.1843]

% A = -A_bar;

% B = [0 0;
%      1 0;
%      0 0; 
%      0 1] + (0.5-rand(4,2));
B_bar = [0 0;
        1 0;
        0 0; 
        0 1];
% B = B_bar;
% B = B_bar + 0.5*(rand(4,2)-0.5);
B = B_bar + 0.5*[ 0.0073   -0.2079;
     -0.1501   -0.1461;
      0.1791   -0.2297;
      0.0391   -0.0779];


C_bar= [1 0 0 0;
    0 0 1 0];
% C = C_bar;
% C = C_bar +0.5*(rand(2,4)-0.5);
C = C_bar + 0.5*[0.2321   -0.2170   -0.1289    0.0539;
            -0.1746   -0.0894   -0.1216   -0.0777];

A_c = kron(eye(N), A);
B_c = kron(eye(N), B);
C_c = kron(eye(N), C);

[~, K, ~] = icare(A_bar,B_bar,3*eye(n),[],[],[],[]); %random stabilizing sol
K = -K;

K_c = kron(eye(N), K);

[~, G, ~] = icare(A_bar',C_bar',3*eye(n),[],[],[],[]); %Should be optimal
G = G';

H = 1; %single agent
re_lambda_list = real(eig(H));
mu = 1.01/min(re_lambda_list);
G = mu*G;
G_c = kron(eye(N), G);

%Design F based on nominal model
F = A_bar+B_bar*K-G*C_bar;


phi_list = 2*pi*rand(1,N); %%initial phases

%derived from sylvester equation
A_f_bar = [A_bar, B_bar*K;
            G*C_bar, kron(eye(N), F)];

sylv_sq_mat = -kron(S', eye(2*n,2*n))+kron(eye(2,2)', A_f_bar);
pz_pl = - kron(eye(2), [C_bar, zeros(2,n)])*inv(sylv_sq_mat)*kron(eye(2,2),[B_bar; zeros(n,2)]); %pd z/pd l based on the nominal model
nominal_pz_pl = pz_pl;

%nominal model IMP
% Pi_bar = zeros(n,2);
% Gamma_bar = zeros(2,2);

% [Pi_bar, Gamma_bar] = IMP_full_solver(A_bar,B_bar,C_bar,zeros(n,2),-eye(2,2), S);
% Sigma_bar = Pi_bar;
%based on nominal model

L_1_bar=reshape(inv(nominal_pz_pl)*Q(:), [2,2]);
L_2_bar= zeros(n,2) ; %Force it to be 0

% L_bar = L_bar - rand(2,2)*3; %debug 

% Pi = zeros(n,2);
% Gamma = zeros(2,2);

% [Pi, Gamma ]= IMP_full_solver(A,B,C,zeros(n,2),-eye(2,2), S);
% L=Gamma - K*Pi

L_1c = kron(eye(N),L_1_bar);
L_2c = kron(eye(N),L_2_bar);

dt = 0.01; %1/dt should be integer
t_max = 4000;
t=0;

X = zeros(N,t_max/dt); % x position
Y = zeros(N,t_max/dt);

% leader_offset = [30; 30];
Omega_c = zeros(2*N,t_max/dt);
X_c =  zeros(n*N,t_max/dt);
Xi_c = zeros(n*N,t_max/dt);

Gramm_c = [];% obtained from [t, t+1]

R_relative = zeros(N,t_max/dt);


X(:,1) =zeros(N,1);
Y(:,1) =zeros(N,1);
X(:,1) =0.3*map_size(1)*(rand(N,1)-0.5);
Y(:,1) =0.3*map_size(2)*(rand(N,1)-0.5);

% X(:,1) =5.3621;
% Y(:,1) = 7.7322;


% X(:,1) =-13.9666;
% Y(:,1) = -1.8377;

% -13.9666
% 0
% -1.8377
% 0

Omega_c(:,1) = v0;
% disp("Initial  positions")
% [X(:,1)'; Y(:,1)']


X_c(1:n:end, 1) = X(:,1); %All other states are set to zero. Xi_c also set to all zero
X_c(2:n:end, 1) = 0;
X_c(3:n:end, 1) = Y(:,1);
X_c(4:n:end, 1) = 0;

disp("Initial conditions")
reshape(X_c(:,1),[n,N])

%actual model
A_closed = [A_c, B_c*K_c, B_c*L_1c;
            kron(H,G*C), kron(eye(N), F), L_2c;
            zeros(2*N, 2*n*N), kron(eye(N), S)];

state_cl = [X_c(:,1); Xi_c(:,1); Omega_c(:,1)]; 


% Prof's idea, optimization like
h_gram = 1; %grammian calculation period
counter = 0; %count the times of grammian calculations
update_L_index = 0;

%For updating the map from l to z
l_list = [];
z_list = []; 

% mapping from l to z based on actual model
A_f_actual = [A, B*K;
            G*C, kron(eye(N), F)];
actual_sylv_sq_mat = -kron(S', eye(2*n,2*n))+kron(eye(2,2)', A_f_actual);
actual_pz_pl = - kron(eye(2), [C, zeros(2,n)])*inv(actual_sylv_sq_mat)*kron(eye(2,2),[B; zeros(n,2)]); 

% L_1_bar=reshape(inv(actual_pz_pl)*Q(:), [2,2]); %debug

if (any(real(eig(A_f_actual)) > -0.05))
    error("actual model is not stable")
end
W_op = diag([100,100,1,1]); % 1/2 (z - q)^T W_op (z-q)
% rt_op = -2*Q(:)';

epsilon_list = []; %store epsilon
epsilon_time = []; %store corresponding time

for k= 2:(t_max/dt)
    t = t+dt;
    state_cl = expm(A_closed*dt)*state_cl;

    X_c(:,k) = state_cl(1:n*N);
    Xi_c(:,k) = state_cl(n*N+1:2*n*N);
    Omega_c(:,k) = state_cl(2*n*N+1:2*n*N+2*N);
        
    for j =1:N
        X(j,k) = X_c(n*(j-1)+1, k); % Assume C = (I_2 0)
        Y(j,k) = X_c(n*(j-1)+3, k);
    end
    if (mod(t,h_gram) < dt) %calculate grammian if (mod(t,1) == 0) does not work, very strange.
        counter = counter+1; 
        k2 = k;
        k1 = k-h_gram/dt;

        temp_gramm = zeros(n*N,2);

       
        for j=(1:N)
            itgl_1 = zeros(2,2); %integral
            itgl_2 = zeros(n,2);
            for l=(k1:k2-1)
                itgl_2 = itgl_2 + X_c(n*(j-1)+1: n*j,l)*Omega_c(2*j-1:2*j,l)'*dt;
                itgl_1 = itgl_1 + Omega_c(2*j-1:2*j,l)* Omega_c(2*j-1:2*j,l)'*dt;
            end
            temp_gramm(n*(j-1)+1:n*j, :) = itgl_2*inv(itgl_1);
        end
        
         
        Gramm_c = [Gramm_c, temp_gramm];

        if (counter > 1 && counter - update_L_index > 1)
            
           
            if (all( abs(Gramm_c(:,end-1:end) - Gramm_c(:,end-3:end-2)) < 0.01) )
                % norm(Gramm_c(:,end-1:end))*0.001
                %debug: get the pi from the actual model
                pi_sigma = -inv(actual_sylv_sq_mat)*kron(eye(2,2), blkdiag(B, eye(n,n)))*reshape([L_1_bar; L_2_bar], [2*(n+2),1]);
                temp_gramm = reshape(pi_sigma, [2*n, 2]);
                temp_gramm = temp_gramm(1:n,:);

                temp_gramm - Gramm_c(:,end-1:end);
                Gramm_c(:,end-1:end) = temp_gramm;
                %debug end

                l_temp = reshape([L_1_bar], [4,1]);
          

                if (update_L_index == 0 || norm(l_temp - l_list(:,end)) > 0.01)
                    disp("pz_pl updated, time: "+ num2str(t))
                    l_list = [l_list, l_temp];
                    z_list = [z_list, reshape(C*Gramm_c(:,end-1:end), [4,1])];
                    
                    temp_mat = kron(l_list', eye(4,4));
                    delta_pz_pl_update = pinv(temp_mat)*(z_list(:) - temp_mat*nominal_pz_pl(:)); % get the solution closest to the nominal one when there are many possible solutions. When the size gets larger, the one get the linear least square sol.
    
                    pz_pl = reshape(nominal_pz_pl(:)+delta_pz_pl_update(:), [4,4]);

                    % % debug
                    % z_temp = actual_pz_pl*l_temp;
                    % z_temp - reshape(C*Gramm_c(:,end-1:end), [4,1])

                    % disp("error of L regression")
                    % norm(z_list(:) - temp_mat*pz_pl(:))
                    % disp("actual error")
                    % norm(z_list(:) - temp_mat*actual_pz_pl(:))
                end
                
                step_length = 1;
                z_1 = reshape(C*Gramm_c(:,end-1:end), [4,1]); % z = vec(C*Pi)
                l_1 = reshape([L_1_bar], [4,1]);
                pe_pz = z_1'*W_op - Q(:)'*W_op;
                pe_pl = pe_pz*pz_pl;
                delta_l = -pe_pl';

                l_2 = l_1 + step_length*delta_l;
                % Pi_2 = reshape(reshape(Pi_1,[2*n,1]) + step_length*delta_pi, [n,2] );
                z_2 = z_1 + pz_pl*(l_2-l_1); %estimated change in z, map: z->l is linear and it is pz_pl (based on nominal model).

                e_1 = 1/2*(z_1-Q(:))'*W_op*(z_1-Q(:));
                e_2 = 1/2*(z_2-Q(:))'*W_op*(z_2-Q(:));

                epsilon_list = [epsilon_list, e_1];
                epsilon_time = [epsilon_time, t];

                % e_2-e_1
                while (e_2-e_1 > 0.3 * step_length* pe_pl * delta_l)
                    step_length = step_length/2;

                    l_2 = l_1 + step_length*delta_l;
                    z_2 = z_1 + pz_pl*(l_2-l_1);

                    e_2 = 1/2*(z_2-Q(:))'*W_op*(z_2-Q(:));
                    % if (step_length < 0.01)
                    %     break;
                    % end
                end

                % while (e_CPi_2-e_CPi_1 > - 1e-2 && e_CPi_2-e_CPi_1 < 0) % make the step bigger to make e_CPi decreases faster
                %     step_length = step_length*1.1
                %     Pi_2 = reshape( reshape(Pi_1,[2*n,1]) - step_length*grad_e, [n,2]);
                %     e_CPi_2 = norm(C*Pi_2-eye(2), 'fro')^2
                % end
    
                %update the L

                d_L = reshape(l_2-l_1, [2,2])
                L_1_bar = L_1_bar + d_L;
                % L_2_bar = L_2_bar + d_L(3:end,:);

                disp("L1 L2 updated, time: "+ num2str(t))
                % L_bar = L_bar + d_L;


                %DEBUG
                % Pi_2
                % temp_sol = linsolve(kron(S', eye(n,n))-kron(eye(2,2)', A+B*K), reshape(B*L_bar, [n*2,1]));

                % Pi_2_actual = reshape(temp_sol, [n,2])
                % e_CPi_2_actual = norm(C*Pi_2_actual-eye(2), 'fro')^2;
                
                %update closed-loop sys
                L_1c = kron(eye(N),L_1_bar);
                % L_2c = kron(eye(N),L_2_bar);
                
                A_closed = [A_c, B_c*K_c, B_c*L_1c;
                            kron(H,G*C), kron(eye(N), F), L_2c;
                            zeros(2*N, 2*n*N), kron(eye(N), S)];

                
    
                update_L_index = counter;

                % W_op = diag([t,t,1,t]);
            end
        end
    end
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end


% rgb_list = rand(N,3);
rgb_list = [255 128 0; 0 153 0; 0 204 204; 204 204 0; 0 51 0; 204 153 0] * 1/255;

E_norm = zeros(N,t_max/dt);

% U_norm = zeros(N,t_max/dt);
for k = 1:t_max/dt
    % u = K_c*X_c(:,k) + L_c*Omega_c(:,k);
    for i = 1:N 
        %Use Omega_1 for all error plottig
        E_norm(i,k) = norm(C*X_c((i-1)*n + 1:  i*n, k) - Q*Omega_c(2*(i-1)+1:2*i ,k));
        % U_norm(i,k) = norm(u((i-1)*2+1: i*2));    
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

% for (i=1:N)
%     disp("U_"+num2str(i)+"^2 * dt/t_max")
%     sum(U_norm(i,:).^2)*dt/t_max
% end

figure(6)
clf(figure(6))
% E = C_c*X_c+Q_c*Omega_; % without using a leader. So the center of circle can be somewhere related to the initial positions of agents
plot(epsilon_time , epsilon_list, '-x','Color', rgb_list(1,:), 'LineWidth', 1);


% figure(1)
% axis equal
% hold on
% clf(figure(1))
% for i=1:N
%     hold on
%     temp_h(i) = plot(X(i,1),Y(i,1),'x', 'MarkerEdgeColor', rgb_list(i,:), 'Color', rgb_list(i,:),'linewidth', 6 );
%     str_array(i) = "Agent"+num2str(i);
%     % disp('s')
% end
% xlabel('x position')
% ylabel('y position')
% legend(temp_h, num2cell(str_array),'AutoUpdate','off');

% sleep_factor = 1; % Prevent the bugs of matlab

% h = [];
% h_agent = []; %moving agent
% for i=1:N
%     h = [h, animatedline( 'Color', rgb_list(i,:))];
%     h_agent = [h_agent, animatedline('Marker', '.', 'MarkerEdgeColor', rgb_list(i,:),'MarkerSize', 20)];
% end

% for k = 1:3*sleep_factor:t_max/dt
%     figure(1)
%     axis equal
%     for i = 1:N
%         addpoints(h(i),X(i,k),Y(i,k));
%         if (k>1)
%             clearpoints(h_agent(i));
%         end
%         addpoints(h_agent(i),X(i,k),Y(i,k));
%     end
%     figure(1)
%     drawnow limitrate nocallbacks
%     pause(dt*sleep_factor);
% end
% drawnow


