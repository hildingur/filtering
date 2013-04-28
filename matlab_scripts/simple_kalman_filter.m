function [ x_hat_minus, x_hat, p, p_minus, K ] = simple_kalman_filter( z, x_init, init_p, A, B, u, Q, R)
% Simple kalman filter, works to estimate state variable x in a the
% following SDE:
% x[k] = A * x[k-1] + B * u[k-1] + w[k-1]
% z[k] = H * x[k] + v[k]
% Where w[k] ~ N(0, Q) and v[k] ~ N(0, R)
% z[k] is directly obervable while x[k] needs to be estimated.
%
%   Detailed explanation goes here
%   z is the measurement vector
%   x_init is the initial guess of the state vector
%   init_p is the initial guess of the posterior estimate covariance
%
%   u is an array of control vectors (We ignore this, it will be set to
%   zero internally, it is only left for compatibility with the paper)
%
%   Q is the co-variance of the state variable x, it will need to be
%   estimated. For our purposes, we know this because we are "cooking" the
%   data.
%
%   R is the co-variance of the state variable z, it will need to be
%   estimated. For our purposes, we know this because we are "cooking" the
%   data.  

    x_hat_minus = zeros(length(z)+1, 1);
    x_hat = zeros(length(z)+1, 1);
    p = zeros(length(z)+1, 1);
    p_minus = zeros(length(z)+1, 1);
    u = zeros(length(z)+1, 1); %Control vector, being ignored
    K = zeros(length(z)+1, 1); %Kalman gain
    H = 1; %The jacobian
    
    x_hat(1) = x_init; %x_hat{^}[k-1] = x_init; The initial guess
    p(1) = init_p;     %P_minus[k-1] = init_p; The initial guess
    
    for k=2:(length(z)+1)
        %
        % The Time update section
        %
        x_hat_minus(k) = A * x_hat(k-1) + B * u(k-1); %eqn 1.9
        p_minus(k) = A * p(k-1) * A + Q; %eqn 1.10, A' = A in 1 dim
        
        %
        % The measurement update section
        %
        K(k) = p_minus(k) * H * (( H * p_minus(k) * H + R) ^ -1); %eqn 1.11
        x_hat(k) = x_hat_minus(k) + K(k) * (z(k-1) - H * x_hat_minus(k)); %eqn 1.12
        p(k) = (1 - K(k) * H) * p_minus(k); %eqn 1.13
    end
    
    
end

