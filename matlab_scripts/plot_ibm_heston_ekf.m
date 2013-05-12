clear;
heston = csvread('../io/ibm_heston_ukf_params.csv', 1);


subplot(2, 2, 1);
omega = heston(:,2);
iterations = heston(:,1);
plot(iterations, omega);
xlim([0, max(iterations)]);
xlabel('# iterations');
ylab = ylabel('\omega', 'rot', 0); %By convention, the second column is omega
last_100_omega = omega(length(omega)-100 : length(omega));
mean = sum(last_100_omega)/length(last_100_omega);
variance = var(last_100_omega);
title(strcat('IBM Heston EKF' ,' last 100 \omega: \mu = ', num2str(mean), ' \sigma^{2} = ', num2str(variance)));


subplot(2, 2, 2);
theta = heston(:,3);
iterations = heston(:,1);
plot(iterations, theta);
xlim([0, max(iterations)]);
xlabel('# iterations');
ylabel('\theta', 'rot', 0); %By convention, the third column is theta
last_100_theta = theta(length(theta)-100 : length(theta));
mean = sum(last_100_theta)/length(last_100_theta);
variance = var(last_100_theta);
title(strcat('IBM Heston EKF' ,' last 100 \theta: \mu = ', num2str(mean), ' \sigma^{2} = ', num2str(variance)));


subplot(2, 2, 3);
xi = heston(:,4);
iterations = heston(:,1);
plot(iterations, xi);
xlim([0, max(iterations)]);
xlabel('# iterations');
ylabel('\xi', 'rot', 0); %By convention, the third column is theta
last_100_xi = xi(length(xi)-100 : length(xi));
mean = sum(last_100_xi)/length(last_100_xi);
variance = var(last_100_xi);
title(strcat('IBM Heston EKF' ,' last 100 \xi: \mu = ', num2str(mean), ' \sigma^{2} = ', num2str(variance)));

subplot(2, 2, 4);
rho = heston(:,5);
iterations = heston(:,1);
plot(iterations, rho);
xlim([0, max(iterations)]);
xlabel('# iterations');
ylabel('\rho', 'rot', 0); %By convention, the third column is theta
last_100_rho = rho(length(rho)-100 : length(rho));
mean = sum(last_100_rho)/length(last_100_rho);
variance = var(last_100_rho);
title(strcat('IBM Heston EKF' ,' last 100 \rho: \mu = ', num2str(mean), ' \sigma^{2} = ', num2str(variance)));

