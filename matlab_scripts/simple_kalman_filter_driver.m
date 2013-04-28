n = 400;  %number of observatoins 
Q = 0.05;    %Variance of the state variable
R = 0.3;    %Variance of the measurement variable


w = normrnd(0, Q, n, 1); %generate a n random variables.

x = 3; %the variable to be measured

x_noisy = x + w; %add noise to x
z = x_noisy + normrnd(0, R, n, 1); % add measurement noise to x_noisy
init_p = 1;
A = 1;
B = 0;
u = 0; %ignored

x_init = 0; %initial guess
% z, x_init, init_p, A, B, u, Q, R
[ x_hat_minus, x_hat, p, p_minus, K ] = simple_kalman_filter(z, x_init, init_p, A, B, u, Q, R);


plot(x_hat, 'Color', 'blue');
title('Simple kalman filter trying to figure out the hidden state variable x = 3');
hold on;
plot(z, 'Color', 'red');

legend('Filtered Signal', 'Original Signal', 'Location','SouthEast');
    
