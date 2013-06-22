%script that summarizes residual results.

subplot(2, 2, 1);
analyze_residuals('./heston_residuals.csv');
title('Heston residuals');

subplot(2, 2, 2);
analyze_residuals('./garch_residuals.csv');
title('GARCH residuals');

subplot(2, 2, 3);
analyze_residuals('./three_two_residuals.csv');
title('3-2 residuals');

subplot(2, 2, 4);
analyze_residuals('./var_p_residuals.csv');
title('p-Model residuals');