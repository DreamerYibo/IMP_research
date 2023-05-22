%%Verify my guess. Can I use IMP to make agents form a circle (maybe a moving circle).

clc; clear;
map_size = [100,100]; % cm
N=3;n=2;
A = zeros(2,2);
B = eye(2,2);
w=1;
l = 3;
S = [0, -w; w, 0];
S_c = kron(eye(N),S);
L_D = [1,-1,0; 0,1,-1; -1,0,1]; %Laplacian of diagraph
K = eye(2); %try this K at first
A_c = kron(eye(N),A)-kron(L_D,B*K);
L_c = [K-S, zeros(2,4); ...
    zeros(2,2), K, zeros(2,2); ...
    K, K-S, zeros(2,2)]; % Try this.
%

dt = 0.01;
t_max = 40;

X = zeros(N,t_max/dt);
Y = zeros(N,t_max/dt);
Omega_ = zeros(2*N,t_max/dt);
X_c =  zeros(2*N,t_max/dt);

R_relative = zeros(N,t_max/dt);

Omega_(:,1) = transpose(l*[1, 0, cos(2/3*pi), sin(2/3*pi), cos(4/3*pi), sin(4/3*pi)]);
X(:,1) = 0.5*map_size(1) + 0.3*map_size(1)*(rand(N,1)-0.5);
Y(:,1) = 0.5*map_size(1) + 0.3*map_size(1)*(rand(N,1)-0.5);
X_c(:,1) = reshape(transpose([X(:,1),Y(:,1)]), [2*N,1]);

for k= 2:(t_max/dt)
    dX_c_dt = A_c*X_c(:,k-1) + L_c*Omega_(:,k-1);
    X_c(:,k) = dX_c_dt*dt+X_c(:,k-1);
    X(:,k) = X_c(1:2:end, k);
    Y(:,k) = X_c(2:2:end, k);

    dOmega_dt = S_c*Omega_(:,k-1);
    Omega_(:,k) =  dOmega_dt*dt + Omega_(:,k-1);

    R(:,k) = ((L_D*X(:,k)).^2 + (L_D*Y(:,k)).^2).^(0.5);
end




rgb_list = rand(N,3);

figure(2)
% plot(transpose(kron(L_D, eye(2))*X_c-kron(eye(N),eye(2))*Omega_))
plot(R')



figure(1)
axis equal
sleep_factor = 1; % Prevent the bugs of matlab

for k = 1:5*sleep_factor:t_max/dt
    figure(1)
    for i = 1:N
        h = animatedline('Marker','.','MarkerEdgeColor', rgb_list(i,:));
        addpoints(h,X(i,k),Y(i,k));
    end
    figure(1)
    drawnow limitrate nocallbacks
    pause(dt*sleep_factor);

end
drawnow

