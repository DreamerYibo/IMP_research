% Unknown A dynamics. B and C are known. Optimization idea. State feedback case. Prof's idea 5.11

% Change A such that the center can be shifted
clc;clear all;
rot_2d = @(t) [cos(t), -sin(t); sin(t), cos(t)];

map_size = [100,100]; % cm
N=1;n=4; %agent 1 original L_bar, agent i L_bar + delta_L_i
w=1;
l = 7;
S = [0, -w; w, 0]; % exp(S_t) is counter clockwise ()
S_bar = zeros(2,2); % Shift the center
exp_St = @(t) [cos(w*t), -sin(w*t); sin(w*t), cos(w*t)];

A_bar = [0 1 0 1;...
     0 0 1 -1;...
     0 0 -1 1;...
     0 0 0 0];

% A = A_bar;

% A = A_bar + rand(n,n)-0.5;
A = A_bar -[ 0.2786   -0.1546    0.1973    0.2723;
            -0.2095   -0.1459   -0.1236    0.2417;
            -0.2604    0.0647    0.2809   -0.1941;
            -0.1072   -0.0933    0.1338   -0.1169];

% A = -A_bar;

% B = [0 0;
%      1 0;
%      0 0; 
%      0 1] + (0.5-rand(4,2));
     B = [0 0;
     1 0;
     0 0; 
     0 1];

C= [1 0 0 0;
    0 0 1 0];

A_c = kron(eye(N), A);
B_c = kron(eye(N), B);
C_c = kron(eye(N), C);

[~, K, ~] = icare(A_bar,B,3*eye(n),[],[],[],[]); %random stabilizing sol
K = -K;

K_c = kron(eye(N), K);



phi_list = 2*pi*rand(1,N); %%initial phases

%nominal model IMP
Pi_bar = zeros(n,2);
Gamma_bar = zeros(2,2);

[Pi_bar, Gamma_bar] = IMP_full_solver(A_bar,B,C,zeros(n,2),-eye(2,2), S);
L_bar=Gamma_bar - K*Pi_bar;

L_bar = L_bar - rand(2,2)*3; %debug 

% Pi = zeros(n,2);
% Gamma = zeros(2,2);

% [Pi, Gamma ]= IMP_full_solver(A,B,C,zeros(n,2),-eye(2,2), S);
% L=Gamma - K*Pi

L_c = kron(eye(N),L_bar);

dt = 0.01; %1/dt should be integer
t_max = 150;
t=0;

X = zeros(N,t_max/dt); % x position
Y = zeros(N,t_max/dt);
% leader_offset = [30; 30];
Omega_c = zeros(2*N,t_max/dt);
X_c =  zeros(n*N,t_max/dt);

Gramm_c = [];% obtained from [t, t+1]

R_relative = zeros(N,t_max/dt);


X(:,1) =zeros(N,1);
Y(:,1) =zeros(N,1);
X(:,1) =0.3*map_size(1)*(rand(N,1)-0.5);
Y(:,1) =0.3*map_size(2)*(rand(N,1)-0.5);
% disp("Initial  positions")
% [X(:,1)'; Y(:,1)']


X_c(1:n:end, 1) = X(:,1); %All other states are set to zero. Xi_c also set to all zero
X_c(2:n:end, 1) = 0;
X_c(3:n:end, 1) = Y(:,1);
X_c(4:n:end, 1) = 0;

disp("Initial conditions")
reshape(X_c(:,1),[n,N])

A_closed = [A_c+ B_c*K_c, B_c*L_c;
           zeros(2*N, n*N), kron(eye(N), S)];

Omega_c(:,1) = reshape([l*cos(phi_list); l*sin(phi_list)], [2*N,1]);

state_cl = [X_c(:,1);Omega_c(:,1)]; 


% Prof's idea, optimization like
h_gram = 1; %grammian calculation period
counter = 0; %count the times of grammian calculations
update_L_index = 0;

%derived from sylvester equation
sylv_sq_mat = kron(S', eye(n,n))-kron(eye(2,2)', A_bar+B*K);
sylv_sq_b = kron(eye(2), B);

% convert to conv optimization without constraint. 1/2* x'Qx + r'x
F_op = linsolve(sylv_sq_mat, sylv_sq_b);
projection_F_op = F_op*inv(F_op'*F_op)*F_op';

Q_op = 2*kron(eye(2,2), C'*C);
rt_op = -2*[[1 0]*C, [0 1]*C];

for k= 2:(t_max/dt)
    t = t+dt;
    state_cl = expm(A_closed*dt)*state_cl;

    X_c(:,k) = state_cl(1:n*N);
    Omega_c(:,k) = state_cl(n*N+1:end);
        
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
            dif_Pi_norm = norm(Gramm_c(:,end-1:end) - Gramm_c(:,end-3:end-2), 'fro');
            if (dif_Pi_norm < 0.01)
                step_length = 1;
                Pi_1 = Gramm_c(:,end-1:end)
                grad_e = Q_op*reshape(Pi_1,[2*n,1]) + rt_op';
                delta_pi = -projection_F_op * grad_e; % direction must be within span(F_op)
                Pi_2 = reshape(reshape(Pi_1,[2*n,1]) + step_length*delta_pi, [n,2] );

                e_CPi_1 = norm(C*Pi_1-eye(2), 'fro')^2
                e_CPi_2 = norm(C*Pi_2-eye(2), 'fro')^2
                while (e_CPi_2-e_CPi_1 > 0.3 * step_length* grad_e' * delta_pi)
                    step_length = step_length/2;
                    Pi_2 = reshape( reshape(Pi_1,[2*n,1]) - step_length*grad_e, [n,2]);
                    e_CPi_2 = norm(C*Pi_2-eye(2), 'fro')^2
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
                d_Pi = Pi_2-Pi_1;
                d_L = inv(-B'*B)*B'*((A_bar+B*K)*d_Pi-d_Pi*S)

                L_bar = L_bar + d_L;

                %DEBUG
                % Pi_2
                % temp_sol = linsolve(kron(S', eye(n,n))-kron(eye(2,2)', A+B*K), reshape(B*L_bar, [n*2,1]));

                % Pi_2_actual = reshape(temp_sol, [n,2])
                % e_CPi_2_actual = norm(C*Pi_2_actual-eye(2), 'fro')^2;
                
                %update closed-loop sys
                L_c = kron(eye(N),L_bar);
                A_closed = [A_c+ B_c*K_c, B_c*L_c;
                zeros(2*N, n*N), kron(eye(N), S)];

                
    
                update_L_index = counter;
            end
        end
    end
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end


% rgb_list = rand(N,3);
rgb_list = [255 128 0; 0 153 0; 0 204 204; 204 204 0; 0 51 0; 204 153 0] * 1/255;

E_norm = zeros(N,t_max/dt);
U_norm = zeros(N,t_max/dt);
for k = 1:t_max/dt
    u = K_c*X_c(:,k) + L_c*Omega_c(:,k);
    for i = 1:N 
        %Use Omega_1 for all error plottig
        E_norm(i,k) = norm(C*X_c((i-1)*n + 1:  i*n, k) - Omega_c(2*(i-1)+1:2*i ,k));
        U_norm(i,k) = norm(u((i-1)*2+1: i*2));
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
 
figure(4)
clf(figure(4))
for (i=1:N)
    temp_h(i) = plot(dt:dt:t_max , U_norm(i,:), 'Color', rgb_list(i,:), 'LineWidth', 1);
    hold on
    str_array(i) = "Agent"+num2str(i);
end
legend(temp_h, num2cell(str_array));
xlabel('time $t \,\, (s)$', 'interpreter', 'latex')
ylabel("$\parallel u_i \,\,\parallel$", 'interpreter', 'latex')

% for (i=1:N)
%     disp("U_"+num2str(i)+"^2 * dt/t_max")
%     sum(U_norm(i,:).^2)*dt/t_max
% end

figure(1)
axis equal
hold on
clf(figure(1))
for i=1:N
    hold on
    temp_h(i) = plot(X(i,1),Y(i,1),'x', 'MarkerEdgeColor', rgb_list(i,:), 'Color', rgb_list(i,:),'linewidth', 6 );
    str_array(i) = "Agent"+num2str(i);
    % disp('s')
end
xlabel('x position')
ylabel('y position')
legend(temp_h, num2cell(str_array),'AutoUpdate','off');

sleep_factor = 1; % Prevent the bugs of matlab

h = [];
h_agent = []; %moving agent
for i=1:N
    h = [h, animatedline( 'Color', rgb_list(i,:))];
    h_agent = [h_agent, animatedline('Marker', '.', 'MarkerEdgeColor', rgb_list(i,:),'MarkerSize', 20)];
end

for k = 1:3*sleep_factor:t_max/dt
    figure(1)
    axis equal
    for i = 1:N
        addpoints(h(i),X(i,k),Y(i,k));
        if (k>1)
            clearpoints(h_agent(i));
        end
        addpoints(h_agent(i),X(i,k),Y(i,k));
    end
    figure(1)
    drawnow limitrate nocallbacks
    pause(dt*sleep_factor);
end
drawnow


