%Have found a strange (A,B) where B has only 1 column.

clc; clear;
N = 1;

A = [0 -1; 2 0];
B = [0; 1];
C = eye(2,2);
Q = eye(2,2);
P = zeros(2,2);

w = 1;
S = [0, -w; w, 0]; % exp(S_t) is counter clockwise ()
exp_St = @(t) [cos(w*t), -sin(w*t); sin(w*t), cos(w*t)];

[~,K,~]=icare(A,B,eye(2,2),1,[],[],[]);
K = -K;

[PI, GAMMA] = mod_IMP_full_solver(A,B,C,P,Q,S);
L = GAMMA-K*PI;

%Add uncertainty on the model
A = A+1e-1*(rand(2,2)-0.5)*A;
B = B+1e-1*(rand(2,2)-0.5)*B;

dt = 0.01;
t_max = 40;

X = zeros(N,t_max/dt); % x position and velocity
Y = zeros(N,t_max/dt);
leader_offset = [30; 30];

Omega_ = zeros(2,t_max/dt);
X_c =  zeros(2*N,t_max/dt);

Omega_(:,1) = [2;0];
Xi_c(:,1) = zeros(4*N,1);
map_size = [100,100]; % cm

X(1) = 0.5*map_size(1) + 0.3*map_size(1)*(rand(N,1)-0.5);
Y(1) = 0.5*map_size(1) + 0.3*map_size(1)*(rand(N,1)-0.5);

X_c(1, 1) = X(1);
X_c(2, 1) = Y(1);

for k= 2:(t_max/dt)
    u = K*X_c(:,k-1)+L*Omega_(:,k-1);
    % u = B_c*K_c*Xi_c(:,k-1);
    dX_c_dt = A*X_c(:,k-1) + B*u;
    X_c(:,k) = dX_c_dt*dt+X_c(:,k-1);
    for j =1:N
        X(j,k) = X_c(2*j-1, k);
        Y(j,k) = X_c(2*j, k);
    end
    Omega_(:,k) = exp_St(k*dt)*Omega_(:,1);
    % R_relative(:,k) = ((L_D*X(1:2:end,k)).^2 + (L_D*Y(1:2:end,k)).^2).^(0.5);
end

rgb_list = rand(N,3);


figure(3)
E = C*X_c+Q*Omega_; % without using a leader. So the center of circle can be somewhere related to the initial positions of agents
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

sleep_factor = 5; % Prevent the bugs of matlab
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

