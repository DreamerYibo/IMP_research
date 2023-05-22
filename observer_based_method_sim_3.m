
% Change A such that the center can be shifted
clc;clear all;
rot_2d = @(t) [cos(t), -sin(t); sin(t), cos(t)];

map_size = [100,100]; % cm
N=4;n=4;
w=1;
l = 7;
S = [0, -w; w, 0]; % exp(S_t) is counter clockwise ()
S_bar = zeros(2,2); % Shift the center
exp_St = @(t) [cos(w*t), -sin(w*t); sin(w*t), cos(w*t)];
% A_D = [0,0,1,0; 0.5,0,0.5,0; 0,0,0,0; 0.5,0,0.5,0];
A_D = [0,0,1,0; 0.5,0,0,0.5; 0,0,0,0.5; 0.5,0,0.5,0];
Delta = diag([0,0,0.5,0]);
% order_list = [3,1,4,2]; % from the leader to the "last" follower
% A_D = zeros(N,N);
H = diag(sum(A_D,2))-A_D+Delta; %Laplacian of diagraph

A = [0 -1 1 2;...
     1 0 3 4;...
     0 0 1 1;...
     0 0 1 0];
B = [0 1; 1 0; 0 1; 0 0];
C= [eye(2), zeros(2,2)];


A_c = kron(eye(N), A);
B_c = kron(eye(N), B);
C_c = kron(eye(N), C);

[~, K, ~] = icare(A,B,3*eye(n),[],[],[],[]); %random stabilizing sol
K = -K;

[~, G, ~] = icare(A',C',B*B',[],[],[],[]); %Should be optimal
G = G';

re_lambda_list = real(eig(H));
mu = 1.01/min(re_lambda_list);
G = mu*G;

K_c = kron(eye(N), K);
G_c = kron(eye(N), G);

D_1 = A+B*K;
D_2 = -G*C;

phi_list = [0,2*pi/3, 4*pi/3, pi]; %%initial phases


Pi = zeros(n,2);
Gamma = zeros(2,2);
Sigma =   [0.8000    0.7900;
        0.1400    0.9600;
        0.4200    0.6600;
        0.9200    0.0400]; %Let Sigma be this such that L_1 L_21 L_22 are all nonzero

[Pi, Gamma] = IMP_full_solver(A,B,C,zeros(n,2),-eye(2,2), S);

% Sigma = Pi; %debug
L_1=Gamma - K*Sigma;
L_21= Sigma*S-D_1*Sigma;
L_22=-(G*C*Pi+D_2*Sigma);
L_1c = kron(eye(N), L_1);
L_2c = kron(eye(N), L_21) + kron(H,L_22);

[Pi, Gamma] = IMP_full_solver(A,B,C,zeros(n,2),-eye(2,2), S_bar);
Sigma = Pi; %Save space in the essay
L_1_bar=Gamma - K*Sigma;
L_21_bar= Sigma*S_bar-D_1*Sigma;
L_22_bar=-(G*C*Pi+D_2*Sigma);
L_1c_bar = kron(eye(N), L_1_bar);
L_2c_bar = kron(eye(N), L_21_bar) + kron(H,L_22_bar);

dt = 0.01;
t_max = 16;

X = zeros(N,t_max/dt); % x position
Y = zeros(N,t_max/dt);
% leader_offset = [30; 30];
Omega_c = zeros(2*N,t_max/dt);
X_c =  zeros(n*N,t_max/dt);
Xi_c = zeros(n*N,t_max/dt);

R_relative = zeros(N,t_max/dt);

Omega_c(:,1) = reshape([l*cos(phi_list); l*sin(phi_list)], [2*N,1]);
Omega_bar = [5;10]; %center
Omega_c_bar = kron(ones(N,1), Omega_bar);
X(:,1) =0.3*map_size(1)*(rand(N,1)-0.5);
Y(:,1) =0.3*map_size(2)*(rand(N,1)-0.5);

% disp("Initial  positions")
% [X(:,1)'; Y(:,1)']


X_c(1:n:end, 1) = X(:,1); %All other states are set to zero. Xi_c also set to all zero
X_c(2:n:end, 1) = Y(:,1);
X_c(3:n:end, 1) = 0;
X_c(4:n:end, 1) = 0;

disp("Initial conditions")
reshape(X_c(:,1),[n,N])

A_closed = [A_c, B_c*K_c, B_c*L_1c, B_c*L_1c_bar;
            kron(H,G*C), kron(eye(N), D_1)+kron(H,D_2), L_2c, L_2c_bar;
            zeros(2*N, 2*n*N), kron(eye(N), S), zeros(2*N, 2*N);
            zeros(2*N, 2*n*N), zeros(2*N, 2*N), kron(eye(N), S_bar)];

state_cl = [X_c(:,1); Xi_c(:,1); Omega_c(:,1); Omega_c_bar(:,1)]; 

for k= 2:(t_max/dt)
    state_cl = expm(A_closed*dt)*state_cl;

    X_c(:,k) = state_cl(1:n*N);
    Xi_c(:,k) = state_cl(n*N+1:2*n*N);
    Omega_c(:,k) = state_cl(2*n*N+1:2*n*N+2*N);
        
    for j =1:N
        X(j,k) = X_c(n*(j-1)+1, k); % Assume C = (I_2 0)
        Y(j,k) = X_c(n*(j-1)+2, k);
    end
    
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end

% rgb_list = rand(N,3);
rgb_list = [255 128 0; 0 153 0; 0 204 204; 204 204 0] * 1/255;

E_norm = zeros(N,t_max/dt);
U_norm = zeros(N,t_max/dt);
for k = 1:t_max/dt
    u = K_c*Xi_c(:,k)+L_1c*Omega_c(:,k)+L_1c_bar*Omega_c_bar;
    for i = 1:N 
        E_norm(i,k) = norm(C*X_c((i-1)*n + 1:  i*n, k) - Omega_c(2*i-1: 2*i ,k)-Omega_bar);
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