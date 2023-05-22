% type 1 controller 
clc;clear;
load("ABK_list");
rot_2d = @(t) [cos(t), -sin(t); sin(t), cos(t)];

map_size = [100,100]; % cm
N=4;n=4;
w=1;
l = 10;
S = [0, -w; w, 0]; % exp(S_t) is counter clockwise ()
exp_St = @(t) [cos(w*t), -sin(w*t); sin(w*t), cos(w*t)];
A_D = [0,0,1,0; 0.5,0,0.5,0; 0,0,0,0; 0.5,0.5,0,0];
% A_D = zeros(N,N);
L_D = eye(N)-A_D; %Laplacian of diagraph

C_1= [eye(2), zeros(2,2)];
C_2 = [zeros(2,2), eye(2)];

C_c = kron(eye(N), C_1);

A_c = [];
B_c = [];
D_K1 = [];
D_K2 = [];
for (i=1:N)
    A_c = blkdiag(A_c, A_list(:,:,i));
    B_c = blkdiag(B_c, B_list(:,:,i));
    D_K1 = blkdiag(D_K1, K_list(:,1:2,i));
    D_K2 = blkdiag(D_K2, K_list(:,3:4,i));
end

K_c = D_K1*(kron(L_D,C_1)) + D_K2*(kron(eye(N),C_2));

Q_c = zeros(2*N,2);
Q_c(1:2,:) = -l*rot_2d(0);
Q_c(3:4,:) = -l*rot_2d(2*pi/3);
Q_c(5:6,:) = -l*rot_2d(4*pi/3);
Q_c(7:8,:) = -l*rot_2d(3*pi/3);

% PI_list = zeros(4,2,N);
% GAMMA_list = zeros(2,2,N);
PI = [];
GAMMA = [];
for (i=1:N)
    [PI_temp, GAMMA_temp] = IMP_full_solver(A_list(:,:,i),B_list(:,:,i),C_1,zeros(4,2),Q_c((2*i-1):(2*i),:), S);
    PI = [PI; PI_temp]; GAMMA = [GAMMA; GAMMA_temp];
end
L_c = GAMMA-K_c*PI;

dt = 0.01;
t_max = 120;

X = zeros(N,t_max/dt); % x position
Y = zeros(N,t_max/dt);
% leader_offset = [30; 30];
Omega_ = zeros(2,t_max/dt);
X_c =  zeros(4*N,t_max/dt);

R_relative = zeros(N,t_max/dt);

Omega_(:,1) = [1;0];
% Omega_(:,1) = transpose(l*[3*1, 0, 2*cos(1/6*pi), 2*sin(1/6*pi), 1/sqrt(3)*cos(1/3*pi), 1/sqrt(3)*sin(1/3*pi)]);
X(:,1) =0.3*map_size(1)*(rand(N,1)-0.5);
Y(:,1) =0.3*map_size(2)*(rand(N,1)-0.5);

X_c(1:4:end, 1) = X(:,1);
X_c(2:4:end, 1) = Y(:,1);
% for k=1:N
%     X_c(4*k-3 : 4*k-2, 1) = X(2*k-1 : 2*k,1) -  leader_offset(1); % Note this
%     X_c(4*k-1 : 4*k, 1) = Y(2*k-1 : 2*k,1) - leader_offset(2);
% end


for k= 2:(t_max/dt)
    u = K_c*X_c(:,k-1)+L_c*Omega_(:,k-1);
    % u = B_c*K_c*Xi_c(:,k-1);
    dX_c_dt = A_c*X_c(:,k-1) + B_c*u;
    X_c(:,k) = dX_c_dt*dt+X_c(:,k-1);
    for j =1:N
        X(j,k) = X_c(4*j-3, k);
        Y(j,k) = X_c(4*j-2, k);
    end
    Omega_(:,k) = exp_St(k*dt)*Omega_(:,1);
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end

rgb_list = rand(N,3);

figure(3)
E = C_c*X_c+Q_c*Omega_; % without using a leader. So the center of circle can be somewhere related to the initial positions of agents
plot(E')
title("Error: C_c*X_c+Q_c*Omega_")

figure(1)
axis equal
hold on
clf(figure(1))
for i=1:N
    hold on
    temp_h(i) = plot(X(i,1),Y(i,1),'x', 'MarkerEdgeColor', rgb_list(i,:), 'linewidth', 10 );
    str_array(i) = "Agent"+num2str(i);
    % disp('s')
end
legend(temp_h, num2cell(str_array),'AutoUpdate','off');

sleep_factor = 10; % Prevent the bugs of matlab
for k = 1:5*sleep_factor:t_max/dt
    figure(1)
    axis equal
    for i = 1:N
        h = animatedline('Marker','.','MarkerEdgeColor', rgb_list(i,:));
        addpoints(h,X(i,k),Y(i,k));
    end
    figure(1)
    drawnow limitrate nocallbacks
    pause(dt*sleep_factor);
end
drawnow