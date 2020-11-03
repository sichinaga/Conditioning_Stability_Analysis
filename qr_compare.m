clear all; close all; clc; 

% Choose initial values of m and n to test 
m0 = 10; 
n0 = 9; 
% Choose number of m and n values to test 
num_m = 10; 
all_m = [1:num_m] * m0;
all_n = [1:num_m] * n0; 
% Choose number of experiments per m, n 
num_exper = 20; 

% Placeholders for data storage 
min_err = zeros(2, num_m);
max_err = zeros(2, num_m);
avg_err = zeros(2, num_m);
all_A = cell(num_m, num_exper); 

% Compute QR error for each mxn matrix 
for i = 1:num_m 
    m = all_m(i);
    n = all_n(i); 
    error1 = zeros(1, num_exper); 
    error2 = zeros(1, num_exper); 
    
    for j = 1:num_exper
        A = randn(m, n); 
        [Q1, R1] = qr(A); 
        [Q2, R2] = qrfactor(A); 
        error1(j) = norm(A - (Q1 * R1)) / norm(A);
        error2(j) = norm(A - (Q2 * R2)) / norm(A);
        all_A{i, j} = A; 
    end 
    
    min_err(1, i) = min(error1);
    max_err(1, i) = max(error1);
    avg_err(1, i) = mean(error1);
    min_err(2, i) = min(error2);
    max_err(2, i) = max(error2);
    avg_err(2, i) = mean(error2);
end 

% Compute QR error for ill-conditioned matrix 
m = max(all_m);
n = max(all_n); 
ill_avg = zeros(1, 2);
ill_A = zeros(num_exper, m, n); 

ill_error1 = zeros(1, num_exper); 
ill_error2 = zeros(1, num_exper);
for j = 1:num_exper
    A = randn(m, n); 
    A(:, end) = A(:, 1); 
    [Q1, R1] = qr(A); 
    [Q2, R2] = qrfactor(A); 
    ill_error1(j) = norm(A - (Q1 * R1)) / norm(A);
    ill_error2(j) = norm(A - (Q2 * R2)) / norm(A);
    ill_A(j, :, :) = A; 
end 
ill_avg(1) = mean(ill_error1);
ill_avg(2) = mean(ill_error2);

% Save results 
save('data.mat', 'all_m', 'all_n', 'num_exper', 'min_err', 'max_err', 'avg_err', 'all_A', 'ill_avg', 'ill_A')