% This code solves (generating relevant plots) the system of ODEs for the first moment vector
% and the distinct entries of the covariance matrix (of a simple gene expression network) for
% 3 cases (positive, zero and negative transcriptional feedback).

%% Initialize
% Time interval

tspan = [0 100];

% Initial conditions

y0 = [ 0 0 0 0 0];

%% Positive Feedback
% Integrating the ODEs

[t,Y] = ode15s(@first_moment_covariance_matrix_pos, tspan, y0);

%%
% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t,Y(:,1));
hold on
plot(t,Y(:,2));
hold off
title('First moments of mRNA and protein copy-numbers')
xlabel('t') 
legend({'y = mRNA copy-numbers','y = Protein copy-numbers'},'Location','southwest')

% Components (distinct) of the covariance matrix
nexttile
plot(t,Y(:,3));
hold on
plot(t,Y(:,4));
hold on
plot(t,Y(:,5));
hold off
title('Covariance matrix entries S_{MM}, S_{MP}, S_{PP}')
xlabel('t') 
legend({'y = S_{MM}','y = S_{MP}','y = S_{PP}'},'Location','southwest')

% Coefficient of variation for mRNA and Protein copy-numbers
CV_mRNA = sqrt(Y(:,3))/Y(:,1);
CV_Protein = sqrt(Y(:,5))/Y(:,2);
nexttile
plot(t,CV_mRNA);
hold on
plot(t,CV_Protein);
hold off
title('Coefficient of variation for mRNA and protein copy-numbers')
xlabel('t') 

%% Zero Feedback
% Integrating the ODEs

[t,Y] = ode15s(@first_moment_covariance_matrix_zero, tspan, y0);

%%
% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t,Y(:,1));
hold on
plot(t,Y(:,2));
hold off
title('First moments of mRNA and protein copy-numbers')
xlabel('t') 
legend({'y = mRNA copy-numbers','y = Protein copy-numbers'},'Location','southwest')

% Components (distinct) of the covariance matrix
nexttile
plot(t,Y(:,3));
hold on
plot(t,Y(:,4));
hold on
plot(t,Y(:,5));
hold off
title('Covariance matrix entries S_{MM}, S_{MP}, S_{PP}')
xlabel('t') 
legend({'y = S_{MM}','y = S_{MP}','y = S_{PP}'},'Location','southwest')

% Coefficient of variation for mRNA and Protein copy-numbers
CV_mRNA = sqrt(Y(:,3))/Y(:,1);
CV_Protein = sqrt(Y(:,5))/Y(:,2);
nexttile
plot(t,CV_mRNA);
hold on
plot(t,CV_Protein);
hold off
title('Coefficient of variation for mRNA and protein copy-numbers')
xlabel('t') 

%% Negative Feedback
% Integrating the ODEs

[t,Y] = ode15s(@first_moment_covariance_matrix_neg, tspan, y0);

%%
% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t,Y(:,1));
hold on
plot(t,Y(:,2));
hold off
title('First moments of mRNA and protein copy-numbers')
xlabel('t') 
legend({'y = mRNA copy-numbers','y = Protein copy-numbers'},'Location','southwest')

% Components (distinct) of the covariance matrix
nexttile
plot(t,Y(:,3));
hold on
plot(t,Y(:,4));
hold on
plot(t,Y(:,5));
hold off
title('Covariance matrix entries S_{MM}, S_{MP}, S_{PP}')
xlabel('t') 
legend({'y = S_{MM}','y = S_{MP}','y = S_{PP}'},'Location','southwest')

% Coefficient of variation for mRNA and Protein copy-numbers
CV_mRNA = sqrt(Y(:,3))/Y(:,1);
CV_Protein = sqrt(Y(:,5))/Y(:,2);
nexttile
plot(t,CV_mRNA);
hold on
plot(t,CV_Protein);
hold off
title('Coefficient of variation for mRNA and protein copy-numbers')
xlabel('t') 
