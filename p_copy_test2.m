% Only make y1 track the given reference. The solutions of IMP equation are not unique in this case
clc;clear all;

map_size = [100,100]; % cm
N=1;n=4; %agent 1 original L_bar, agent i L_bar + delta_L_i
w=1;
l = 7;
S = [0, 1; 0, 0]; % exp(S_t) is counter clockwise ()
Q = [1 -1]; % y1 ---> Qv. y ---> [t; t]
v0 = [0 1];


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
A = A_bar +  [0.0573   -0.2749         0    0.1002;
            0.3933   -0.3173         0    0.2401;
            -0.4619   -0.3687         0   -0.3001;
            0.4954    0.4156         0    0.0686;];


% A = -A_bar;

% B = [0 0;
%      1 0;
%      0 0; 
%      0 1] + (0.5-rand(4,2));
B_bar = [0 0;
        1 0;
        0 0; 
        0 1];
B = B_bar;
% B = B_bar + [ 0.0073   -0.2079;
%      -0.1501   -0.1461;
%       0.1791   -0.2297;
%       0.0391   -0.0779];

C_bar= [1 0 1 0;];

C = C_bar;
% C = C_bar + [0.2321   -0.2170   -0.1289    0.0539;];


A_c = kron(eye(N), A);
B_c = kron(eye(N), B);
C_c = kron(eye(N), C);

% p-copy
G_1 = kron(eye(1), [0 1; 0 0]);
G_2 = kron(eye(1), [0; 1]);

[~, K, ~] = icare([A_bar, zeros(n,2); G_2*C_bar, G_1], [B_bar; zeros(2,2)] ,3*eye(n+2),[],[],[],[]); %random stabilizing sol
K = -K;

%%some other tests: start

n = size(A,1); % size of x
m = size(B,2); % size of u
r = size(S,2); %size of omega
p = size(C,1); %size of e

A1 = [A, B; C, zeros(p,m)];
    temp = zeros(n+p,n+m);
    temp(1:n, 1:n) = eye(n,n);
    A2 = temp;

    temp2 = kron(transpose(eye(r,r)), A1) -  kron(transpose(S), A2);
    temp2 = kron(transpose(eye(r,r)), A1) -  kron(transpose(S), A2);
    null_temp2 = null(temp2)
    null_pi_gamma = reshape(null_temp2, [4+2,2,size(null_temp2, 2)])

%%some other tests: end




dt = 0.01;
t_max = 300;

X = zeros(N,t_max/dt); % x position
Y = zeros(N,t_max/dt);
% leader_offset = [30; 30];
Omega_c = zeros(2*N,t_max/dt);
X_c =  zeros(n*N,t_max/dt);
Xi_c = zeros(2*N,t_max/dt);

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

A_closed = [[A, zeros(n,2); G_2*C, G_1] + [B; zeros(2,2)]*K, [zeros(n,2); G_2*(-Q)];
            zeros(2, n+2), kron(eye(N), S);];

state_cl = [X_c(:,1); Xi_c(:,1); Omega_c(:,1)]; 

for k= 2:(t_max/dt)
    state_cl = expm(A_closed*dt)*state_cl;

    X_c(:,k) = state_cl(1:n);
    Xi_c(:,k) = state_cl(n+1:n+2);
    Omega_c(:,k) = state_cl(n+3:end);
        
    for j =1:N
        X(j,k) = X_c(n*(j-1)+1, k); % Assume C = (I_2 0)
        Y(j,k) = X_c(n*(j-1)+3, k);
    end
    
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end

% rgb_list = rand(N,3);
rgb_list = [255 128 0; 0 153 0; 0 204 204; 204 204 0] * 1/255;

E_norm = zeros(N,t_max/dt);
U_norm = zeros(N,t_max/dt);
for k = 1:t_max/dt
    u = K*[X_c(:,k);Xi_c(:,k)];
    for i = 1:N 
        E_norm(i,k) = norm(C*X_c((i-1)*n + 1:  i*n, k) - Q*Omega_c(2*i-1: 2*i ,k));
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

% %DEBUG
% U_norm_sdy = zeros(N,t_max/dt); %steady state u_norm
% error_sim = zeros(1,t_max/dt);
% figure(7)
% clf(figure(7))
% for k = 1:t_max/dt
%     u = K_c*Xi_c(:,k)+L_1c*Omega_c(:,k)+L_1c_bar*Omega_c_bar;
%     u_sdy = kron(eye(N),Gamma_1)*Omega_c(:,k)+kron(eye(N),Gamma_2)*Omega_c_bar;
%     %debug
%     error_sim(k) = norm(u-u_sdy);
%     for i = 1:N 
%         U_norm_sdy(i,k) = norm(u_sdy((i-1)*2+1: i*2));
%     end
% end
% for (i=1:N)
%     temp_h(i) = plot(dt:dt:t_max , U_norm_sdy(i,:), 'Color', rgb_list(i,:), 'LineWidth', 1);
%     hold on
%     str_array(i) = "Agent"+num2str(i);
% end

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

