function [residuals, csvfile] = analyze_residuals(filename)
%filename = './heston_residuals.csv';

%filename = 'C:/Users/hildingur/Downloads/StdResiduals.csv';

%Reading in the csv file.
csvfile = csvread(filename, 1, 4);
residuals = csvfile(:,1);

mpe = mean(residuals);
rmse = std(residuals);

disp(strcat('MPE = ', num2str(mpe)));
disp(strcat('RMSE = ', num2str(rmse)));

%normalizing the residuals.
residuals = (residuals - mean(residuals)) / std(residuals);

% check to see if the data passes the Chi-squared gof test
% for normalcy.
[h,p,st] = chi2gof(residuals, 'nbins', 20, 'cdf', @normcdf);
disp(strcat('chi square statistics are h = ', num2str(h) ...
    , '  p = ', num2str(p),...
    ' st.df = ', num2str(st.df), ...
    ' st.chi2stat = ', num2str(st.chi2stat)));

if h==0
    disp('Residuals are normal');
else
    disp(strcat('Residuals are NOT normal within the 95% confidence interval'));
end

% check to see if the residuals pass the Ljung-Box test for
% auto-correlation
h = lbqtest(residuals);
if h==0
    disp('Residuals are NOT auto-correlated');
else
    disp(strcat('Residuals are auto-correlated'));
end



%plot the normalized histogram of the resiiduals vs the normal pdf
histnorm(residuals, 100);
hold on;
x = -6:.1:6;
plot(x, normpdf(x, 0, 1), 'r-');
axis([-10, 10, 0, .6]);
set(gca,'XTick',-10:2:10);
hold off;

ylabel('Frequency');
xlabel('Errors');