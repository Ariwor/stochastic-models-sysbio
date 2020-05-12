% Analysis of the Schlogl model with rate constants that depend on a volume parameter V, which is defined by the user.
% Two different moment closure methods are employed, Mass Fluctuation Kinetics(MFK) and Separable Derivative Matching(SDM) and
% are compared with the results obtained from the Stochastic Simulation Algorithm (SSA).

%%
% Time interval
tspan = [0 50];

% Initial conditions

y0 = [ 0 0 0 0 0];
y0_perturbed = y0 + 10^(-1); %Perturb initial values for SDM integration

%% Mass Fluctuation Kinetics(MFK) 

% Integrating the ODEs

V = [2, 10];
[t_MFK_2,Y_MFK_2] = ode15s(@first_moment_covariance_MFK_Schlogl, tspan, y0, [], V(1));
[t_MFK_10,Y_MFK_10] = ode15s(@first_moment_covariance_MFK_Schlogl, tspan, y0, [], V(2));

%% Separable Derivative Matching(SDM)

% Integrating the ODEs

[t_SDM_2,Y_SDM_2] = ode15s(@first_moment_covariance_SDM_Schlogl, tspan, y0_perturbed, [], V(1));
[t_SDM_10,Y_SDM_10] = ode15s(@first_moment_covariance_SDM_Schlogl, tspan, y0_perturbed, [], V(2));

%% SSA

T_final = 50;
V_2 = 2;
V_10 = 10;
P = 10000; % Number of paths
mRNA_SSA_vol_2=zeros(P,length(t_MFK_2));
Protein_SSA_vol_2=zeros(P,length(t_MFK_2));
mRNA_SSA_vol_10=zeros(P,length(t_MFK_10));
Protein_SSA_vol_10=zeros(P,length(t_MFK_10));

% Volume = 2

for k=1:P
[X_results_vol_2, jumps_vol_2] = SSA_Schlogl_vol(T_final, V_2);
    for index=1:length(jumps_vol_2)
        first_jump_index=find(t_MFK_2>jumps_vol_2(index),1,'first');
        mRNA_SSA_vol_2(k,first_jump_index:end)=X_results_vol_2(index, 1);
        Protein_SSA_vol_2(k,first_jump_index:end)=X_results_vol_2(index, 2);
    end
end

% Volume = 10

for k=1:P
[X_results_vol_10, jumps_vol_10] = SSA_Schlogl_vol(T_final, V_10);
    for index=1:length(jumps_vol_10)
        first_jump_index=find(t_MFK_10>jumps_vol_10(index),1,'first');
        mRNA_SSA_vol_10(k,first_jump_index:end)=X_results_vol_10(index, 1);
        Protein_SSA_vol_10(k,first_jump_index:end)=X_results_vol_10(index, 2);
    end
end

%% Moments from SSA

% Volume = 2

y_SSA_vol_2=zeros(5,length(t_MFK_2));
for index_SSA=1:length(t_MFK_2)
    y_SSA_vol_2(1, index_SSA)=mean(mRNA_SSA_vol_2(:,index_SSA));
    y_SSA_vol_2(2,index_SSA)=mean(Protein_SSA_vol_2(:,index_SSA));
    Sigma_cov=cov(mRNA_SSA_vol_2(:,index_SSA),Protein_SSA_vol_2(:,index_SSA));
    y_SSA_vol_2(3,index_SSA)=Sigma_cov(1,1);
    y_SSA_vol_2(4,index_SSA)=Sigma_cov(1,2);
    y_SSA_vol_2(5,index_SSA)=Sigma_cov(2,2);
end
CV_mRNA_SSA_vol_2 = sqrt(y_SSA_vol_2(3,:))./y_SSA_vol_2(1,:);
CV_Protein_SSA_vol_2 = sqrt(y_SSA_vol_2(5,:))./y_SSA_vol_2(2,:);


% Volume = 10

y_SSA_vol_10=zeros(5,length(t_MFK_10));
for index_SSA=1:length(t_MFK_10)
    y_SSA_vol_10(1, index_SSA)=mean(mRNA_SSA_vol_10(:,index_SSA));
    y_SSA_vol_10(2,index_SSA)=mean(Protein_SSA_vol_10(:,index_SSA));
    Sigma_cov=cov(mRNA_SSA_vol_10(:,index_SSA),Protein_SSA_vol_10(:,index_SSA));
    y_SSA_vol_10(3,index_SSA)=Sigma_cov(1,1);
    y_SSA_vol_10(4,index_SSA)=Sigma_cov(1,2);
    y_SSA_vol_10(5,index_SSA)=Sigma_cov(2,2);
end
CV_mRNA_SSA_vol_10 = sqrt(y_SSA_vol_10(3,:))./y_SSA_vol_10(1,:);
CV_Protein_SSA_vol_10 = sqrt(y_SSA_vol_10(5,:))./y_SSA_vol_10(2,:);

%% Plotting (Volume = 2)
%% Plotting MFK - SSA

% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t_MFK_2,Y_MFK_2(:,1));
hold on
plot(t_MFK_2,y_SSA_vol_2(1,:));
plot(t_MFK_2,Y_MFK_2(:,2));
plot(t_MFK_2,y_SSA_vol_2(2,:));
hold off
title('First moments of X and Y from MFK (volume = 2)')
xlabel('t') 
legend({'y = E_{X} from MFK','y = E_{X} from SSA','y = E_{Y} from MFK','y = E_{Y} from SSA'},'Location','northeast')

% Components (distinct) of the covariance matrix
nexttile
plot(t_MFK_2,Y_MFK_2(:,3));
hold on
plot(t_MFK_2,y_SSA_vol_2(3,:));
plot(t_MFK_2,Y_MFK_2(:,4));
plot(t_MFK_2,y_SSA_vol_2(4,:));
plot(t_MFK_2,Y_MFK_2(:,5));
plot(t_MFK_2,y_SSA_vol_2(5,:));
hold off
title('Covariance matrix entries S_{XX}, S_{XY}, S_{YY} from MFK (volume = 2)')
xlabel('t') 
legend({'y = S_{XX} from MFK','y = S_{XX} from SSA','y = S_{XY} from MFK','y = S_{XY} from SSA', 'y = S_{YY} from MFK','y = S_{YY} from SSA'},'Location','northeast')

% Coefficient of variation for X and Y
CV_mRNA_MFK_vol_2 = sqrt(Y_MFK_2(:,3))./Y_MFK_2(:,1);
CV_Protein_MFK_vol_2 = sqrt(Y_MFK_2(:,5))./Y_MFK_2(:,2);
nexttile
plot(t_MFK_2,CV_mRNA_MFK_vol_2);
hold on
plot(t_MFK_2,CV_mRNA_SSA_vol_2);
plot(t_MFK_2,CV_Protein_MFK_vol_2);
plot(t_MFK_2,CV_Protein_SSA_vol_2);
hold off
title('Coefficient of variation for X and Y from MFK (volume = 2)')
xlabel('t')
legend({'y = CV_{X} from MFK','y = CV_{X} from SSA','y = CV_{Y} from MFK','y = CV_{Y} from SSA'},'Location','northeast')

%% Plotting SDM - SSA

% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t_SDM_2,Y_SDM_2(:,1));
hold on
plot(t_MFK_2,y_SSA_vol_2(1,:));
plot(t_SDM_2,Y_SDM_2(:,2));
plot(t_MFK_2,y_SSA_vol_2(2,:));
hold off
title('First moments of X and Y from SDM (volume = 2)')
xlabel('t') 
legend({'y = E_{X} from SDM','y = E_{X} from SSA','y = E_{Y} from SDM','y = E_{Y} from SSA'},'Location','northeast')
%
% Components (distinct) of the covariance matrix
nexttile
plot(t_SDM_2,Y_SDM_2(:,3));
hold on
plot(t_MFK_2,y_SSA_vol_2(3,:));
plot(t_SDM_2,Y_SDM_2(:,4));
plot(t_MFK_2,y_SSA_vol_2(4,:));
plot(t_SDM_2,Y_SDM_2(:,5));
plot(t_MFK_2,y_SSA_vol_2(5,:));
hold off
title('Covariance matrix entries S_{XX}, S_{XY}, S_{YY} from SDM (volume = 2)')
xlabel('t') 
legend({'y = S_{XX} from SDM','y = S_{XX} from SSA','y = S_{XY} from SDM','y = S_{XY} from SSA', 'y = S_{YY} from SDM','y = S_{YY} from SSA'},'Location','northeast')

% Coefficient of variation for X and Y
CV_mRNA_SDM_vol_2 = sqrt(Y_SDM_2(:,3))./Y_SDM_2(:,1);
CV_Protein_SDM_vol_2 = sqrt(Y_SDM_2(:,5))./Y_SDM_2(:,2);
nexttile
plot(t_SDM_2,CV_mRNA_SDM_vol_2);
hold on
plot(t_MFK_2,CV_mRNA_SSA_vol_2);
plot(t_SDM_2,CV_Protein_SDM_vol_2);
plot(t_MFK_2,CV_Protein_SSA_vol_2);
hold off
title('Coefficient of variation for X and Y from SDM (volume = 2)')
xlabel('t')
legend({'y = CV_{X} from SDM','y = CV_{X} from SSA','y = CV_{Y} from SDM','y = CV_{Y} from SSA'},'Location','northeast')

%% Plotting (Volume = 10)
%% Plotting MFK - SSA

% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t_MFK_10,Y_MFK_10(:,1));
hold on
plot(t_MFK_10,y_SSA_vol_10(1,:));
plot(t_MFK_10,Y_MFK_10(:,2));
plot(t_MFK_10,y_SSA_vol_10(2,:));
hold off
title('First moments of X and Y from MFK (volume = 10)')
xlabel('t') 
legend({'y = E_{X} from MFK','y = E_{X} from SSA','y = E_{Y} from MFK','y = E_{Y} from SSA'},'Location','northeast')

% Components (distinct) of the covariance matrix
nexttile
plot(t_MFK_10,Y_MFK_10(:,3));
hold on
plot(t_MFK_10,y_SSA_vol_10(3,:));
plot(t_MFK_10,Y_MFK_10(:,4));
plot(t_MFK_10,y_SSA_vol_10(4,:));
plot(t_MFK_10,Y_MFK_10(:,5));
plot(t_MFK_10,y_SSA_vol_10(5,:));
hold off
title('Covariance matrix entries S_{XX}, S_{XY}, S_{YY} from MFK (volume = 10)')
xlabel('t') 
legend({'y = S_{XX} from MFK','y = S_{XX} from SSA','y = S_{XY} from MFK','y = S_{XY} from SSA', 'y = S_{YY} from MFK','y = S_{YY} from SSA'},'Location','northeast')

% Coefficient of variation for X and Y
CV_mRNA_MFK_vol_10 = sqrt(Y_MFK_10(:,3))./Y_MFK_10(:,1);
CV_Protein_MFK_vol_10 = sqrt(Y_MFK_10(:,5))./Y_MFK_10(:,2);
nexttile
plot(t_MFK_10,CV_mRNA_MFK_vol_10);
hold on
plot(t_MFK_10,CV_mRNA_SSA_vol_10);
plot(t_MFK_10,CV_Protein_MFK_vol_10);
plot(t_MFK_10,CV_Protein_SSA_vol_10);
hold off
title('Coefficient of variation for X and Y from MFK (volume = 10)')
xlabel('t')
legend({'y = CV_{X} from MFK','y = CV_{X} from SSA','y = CV_{Y} from MFK','y = CV_{Y} from SSA'},'Location','northeast')

%% Plotting SDM - SSA

% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t_SDM_10,Y_SDM_10(:,1));
hold on
plot(t_MFK_10,y_SSA_vol_10(1,:));
plot(t_SDM_10,Y_SDM_10(:,2));
plot(t_MFK_10,y_SSA_vol_10(2,:));
hold off
title('First moments of X and Y from SDM (volume = 10)')
xlabel('t') 
legend({'y = E_{X} from SDM','y = E_{X} from SSA','y = E_{Y} from SDM','y = E_{Y} from SSA'},'Location','northeast')

% Components (distinct) of the covariance matrix
nexttile
plot(t_SDM_10,Y_SDM_10(:,3));
hold on
plot(t_MFK_10,y_SSA_vol_10(3,:));
plot(t_SDM_10,Y_SDM_10(:,4));
plot(t_MFK_10,y_SSA_vol_10(4,:));
plot(t_SDM_10,Y_SDM_10(:,5));
plot(t_MFK_10,y_SSA_vol_10(5,:));
hold off
title('Covariance matrix entries S_{XX}, S_{XY}, S_{YY} from SDM (volume = 10)')
xlabel('t') 
legend({'y = S_{XX} from SDM','y = S_{XX} from SSA','y = S_{XY} from SDM','y = S_{XY} from SSA', 'y = S_{YY} from SDM','y = S_{YY} from SSA'},'Location','northeast')

% Coefficient of variation for X and Y
CV_mRNA_SDM_vol_10 = sqrt(Y_SDM_10(:,3))./Y_SDM_10(:,1);
CV_Protein_SDM_vol_10 = sqrt(Y_SDM_10(:,5))./Y_SDM_10(:,2);
nexttile
plot(t_SDM_10,CV_mRNA_SDM_vol_10);
hold on
plot(t_MFK_10,CV_mRNA_SSA_vol_10);
plot(t_SDM_10,CV_Protein_SDM_vol_10);
plot(t_MFK_10,CV_Protein_SSA_vol_10);
hold off
title('Coefficient of variation for X and Y from SDM (volume = 10)')
xlabel('t')
legend({'y = CV_{X} from SDM','y = CV_{X} from SSA','y = CV_{Y} from SDM','y = CV_{Y} from SSA'},'Location','northeast')