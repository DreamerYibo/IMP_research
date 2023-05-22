
% Observer Based method. Able to know neighbor's observer states.
clc;clear;
rot_2d = @(t) [cos(t), -sin(t); sin(t), cos(t)];

map_size = [100,100]; % cm
N=4;n=3;
w=1;
l = 7;
S = [0, -w; w, 0]; % exp(S_t) is counter clockwise ()
exp_St = @(t) [cos(w*t), -sin(w*t); sin(w*t), cos(w*t)];
% A_D = [0,0,1,0; 0.5,0,0.5,0; 0,0,0,0; 0.5,0,0.5,0];
A_D = [0,0,1,0; 0.5,0,0.5,0; 0,0,0,0.5; 0.5,0,0.5,0];
order_list = [3,1,4,2]; % from the leader to the "last" follower
% A_D = zeros(N,N);
L_D = eye(N)-A_D; %Laplacian of diagraph

A = [0 -1 1; 1 0 0; 2 0 1];
B = [1 0; 0 0; 0 1];
C= [eye(2), zeros(2,1)];


A_c = kron(eye(N), A);
B_c = kron(eye(N), B);
C_c = kron(eye(N), C);

[~, K, ~] = icare(A,B,C'*C,[],[],[],[]); %random stabilizing sol
K = -K;

[~, G, ~] = icare(A',C',B*B',[],[],[],[]); %Should be optimal
G = G';

re_lambda_list = real(eig(L_D));
mu = 1.01/min(re_lambda_list);
G = mu*G;

K_c = kron(eye(N), K);
G_c = kron(eye(N), G);


Q_c = zeros(2*N,2);
Q_c(1:2,:) = -l*rot_2d(0);
Q_c(3:4,:) = -l*rot_2d(2*pi/3);
Q_c(5:6,:) = -l*rot_2d(4*pi/3);
Q_c(7:8,:) = -l*rot_2d(3*pi/3);

PI = [];
GAMMA = [];

for (i=1:N)
    [PI_temp, GAMMA_temp] = IMP_full_solver(A,B,C,zeros(n,2),Q_c((2*i-1):(2*i),:), S);
    PI = [PI; PI_temp]; GAMMA = [GAMMA; GAMMA_temp];
end
L_c = GAMMA-K_c*PI;

dt = 0.01;
t_max = 20;

X = zeros(N,t_max/dt); % x position
Y = zeros(N,t_max/dt);
% leader_offset = [30; 30];
Omega_ = zeros(2,t_max/dt);
X_c =  zeros(n*N,t_max/dt);
Xi_c = zeros(n*N,t_max/dt);

R_relative = zeros(N,t_max/dt);

Omega_(:,1) = [1;0];
X(:,1) =0.3*map_size(1)*(rand(N,1)-0.5);
Y(:,1) =0.3*map_size(2)*(rand(N,1)-0.5);

X_c(1:n:end, 1) = X(:,1);
X_c(2:n:end, 1) = Y(:,1);


for k= 2:(t_max/dt)
    u = K_c*Xi_c(:,k-1)+L_c*Omega_(:,k-1);
    % u = B_c*K_c*Xi_c(:,k-1);
    dX_c_dt = A_c*X_c(:,k-1) + B_c*u;
    dXi_c_dt = (A_c-G_c*C_c*kron(L_D, eye(n)))*Xi_c(:,k-1)+G_c*C_c*(kron(L_D, eye(n)))*X_c(:,k-1)+B_c*u;

    X_c(:,k) = dX_c_dt*dt+X_c(:,k-1);
    Xi_c(:,k) = dXi_c_dt*dt+Xi_c(:,k-1);
        
    for j =1:N
        X(j,k) = X_c(n*(j-1)+1, k);
        Y(j,k) = X_c(n*(j-1)+2, k);
    end
    Omega_(:,k) = exp_St(k*dt)*Omega_(:,1);
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end

rgb_list = rand(N,3);

E_norm = zeros(N,t_max/dt);
for k = 1:t_max/dt
    for i = 1:N 
        E_norm(i,k) = norm(C*X_c((i-1)*n + 1:  i*n, k) + Q_c((i-1)*2+1:i*2, :)*Omega_(:,k));
    end
end
figure(3)
clf(figure(3))
% E = C_c*X_c+Q_c*Omega_; % without using a leader. So the center of circle can be somewhere related to the initial positions of agents
for (i=1:N)
    temp_h(i) = plot(dt:dt:t_max , E_norm(i,:), 'Color', rgb_list(i,:), 'LineWidth', 1);
    hold on
    str_array(i) = "Agent"+num2str(i);
end
legend(temp_h, num2cell(str_array));
xlabel('time $t \,\, (s)$', 'interpreter', 'latex')
title("$\parallel Cx_i-Q_i v \,\,\parallel$", 'interpreter', 'latex')

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