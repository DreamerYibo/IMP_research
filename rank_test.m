clc; clear;


n = 12;
N = n*10;
matrix = [];
v = [];
for i=1:N 
    for j=1:n 
        phi = 2*pi*rand;
        v = [v; cos(phi); sin(phi)];
    end
    matrix = [matrix, v];
    v = [];
end
