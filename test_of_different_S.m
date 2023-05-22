clc; clear;

load ABK_list.mat;

A = A_list(:,:,4);
B = B_list(:,:,4);
Mat_solva = @(t) [A+1i*t*eye(4,4), B; C_1, zeros(2,2)];
N= 1e5;
sigma_min_list = zeros(1,N);
det_list = zeros(1,N);
k=1;
for t = linspace(0,1000, N)
    [~, Sigma, ~] = svd(Mat_solva(t));
    sigma_min_list(k) = Sigma(end,end);
    det_list(k) = det(Mat_solva(t));
    k = k+1;
end

figure(1)
hold on
title("|det_list|")
plot(abs(det_list), 'o');
hold off

figure(2)
hold on
title("sigma_min_list")
plot(sigma_min_list);
hold off