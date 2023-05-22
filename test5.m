%% Try making follower use only the relative position. (With leader)
clc; clear;

rot_2d = @(t) [cos(t), -sin(t); sin(t), cos(t)];

map_size = [100,100]; % cm
N=3;n=4;
A = zeros(4,4);
A(1,2)=1;A(3,4)=1;
B = zeros(4,2);
B(2,1)=1;B(4,2)=1;
w=1;
l = 7;
S = [0, -w; w, 0]; % exp(S_t) is counter clockwise ()
exp_St = @(t) [cos(w*t), -sin(w*t); sin(w*t), cos(w*t)];
L_D = [1,-1,0; 0,1,-1; -1,0,1]; %Laplacian of diagraph

k_pos=-1;
k_v = -2;
K1=zeros(2,4);
K1(1,1) = k_pos; K1(1,2)=k_v; K1(2,3) = k_pos; K1(2,4)=k_v;
K2=zeros(2,4);
K2(1,1) = k_pos; K2(2,3)=k_pos;
K3=K1;

A_c = kron(eye(N),A);
B_c = kron(eye(N),B);
% K_c = [K1,-K2,zeros(2,4); zeros(2,4), K1,-K2; -K2,zeros(2,4), K1]; % Note that there is no leader 
K_c = [K1,-K2,zeros(2,4); zeros(2,4), K1,-K2; zeros(2,8), K3]; % Do this to make sure A_c - B_c*K_c is asymptotic stable

C = zeros(2,4);
C(1,1)=1;C(2,3)=1;
C_c = kron(eye(3), C);

% internal states
G = -transpose(K1); %try this 3
Gc_1 = kron(eye(3), G*C);
Gc_2 = kron([1,-1,0;0,1,-1; 0,0,1], G*C);

Q_c = zeros(6,2);
Q_c(1:2,:) = -l*rot_2d(0);
Q_c(3:4,:) = -l*rot_2d(2*pi/3);
Q_c(5:6,:) = -l*rot_2d(4*pi/3);


[PI, GAMMA] = IMP_full_solver(A_c,B_c,C_c,zeros(12,2),Q_c,S);
L_c = GAMMA-K_c*PI;

dt = 0.01;
t_max = 120;

X = zeros(2*N,t_max/dt); % x position and velocity
Y = zeros(2*N,t_max/dt);
leader_offset = [30; 30];
Xi_c = zeros(4*N,t_max/dt);
Omega_ = zeros(2,t_max/dt);
X_c =  zeros(4*N,t_max/dt);

R_relative = zeros(N,t_max/dt);

Omega_(:,1) = [1;0];
Xi_c(:,1) = zeros(4*N,1);
% Omega_(:,1) = transpose(l*[3*1, 0, 2*cos(1/6*pi), 2*sin(1/6*pi), 1/sqrt(3)*cos(1/3*pi), 1/sqrt(3)*sin(1/3*pi)]);
X(1:2:end,1) = 0.5*map_size(1) + 0.3*map_size(1)*(rand(N,1)-0.5);
Y(1:2:end,1) = 0.5*map_size(1) + 0.3*map_size(1)*(rand(N,1)-0.5);

for k=1:N
    X_c(4*k-3 : 4*k-2, 1) = X(2*k-1 : 2*k,1) -  leader_offset(1); % Note this
    X_c(4*k-1 : 4*k, 1) = Y(2*k-1 : 2*k,1) - leader_offset(2);
end


for k= 2:(t_max/dt)
    u = B_c*(K_c-kron(eye(3), K1-K2))*X_c(:,k-1) + B_c*(kron(eye(3), K1-K2))*Xi_c(:,k-1) + B_c*L_c*Omega_(:,k-1);
    % u = B_c*K_c*Xi_c(:,k-1);
    dX_c_dt = A_c*X_c(:,k-1) + u;
    dXi_c_dt = (A_c-Gc_1)*Xi_c(:,k-1)+Gc_2*X_c(:,k-1)+ u - [G*Q_c(3:4,:)*Omega_(:,k-1); G*Q_c(5:6,:)*Omega_(:,k-1); zeros(4,1)]; % try this
    X_c(:,k) = dX_c_dt*dt+X_c(:,k-1);
    Xi_c(:,k) = dXi_c_dt*dt+Xi_c(:,k-1);
    for j =1:N
        X(2*j-1 : 2*j,k) = X_c(4*j-3 : 4*j-2, k) + leader_offset(1); % Note this
        Y(2*j-1 : 2*j,k) = X_c(4*j-1 : 4*j, k) + leader_offset(2);
    end
    % dOmega_dt = S_c*Omega_(:,k-1);
    % Omega_(:,k) =  dOmega_dt*dt + Omega_(:,k-1);
    Omega_(:,k) = exp_St(k*dt)*Omega_(:,1);
    R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end




rgb_list = rand(N,3);

figure(2)
% plot(transpose(kron(L_D, eye(2))*X_c-kron(eye(N),eye(2))*Omega_))
plot(R_relative')
title("Relative distance between agents")

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
    temp_h(i) = plot(X(2*i-1,1),Y(2*i-1,1),'x', 'MarkerEdgeColor', rgb_list(i,:), 'linewidth', 10 );
    str_array(i) = "Agent"+num2str(i);
    % disp('s')
end
legend(temp_h, num2cell(str_array),'AutoUpdate','off');

sleep_factor = 1; % Prevent the bugs of matlab
for k = 1:5*sleep_factor:t_max/dt
    figure(1)
    axis equal
    for i = 1:N
        h = animatedline('Marker','.','MarkerEdgeColor', rgb_list(i,:));
        addpoints(h,X(2*i-1,k),Y(2*i-1,k));
    end
    figure(1)
    drawnow limitrate nocallbacks
    pause(dt*sleep_factor);

end
drawnow

