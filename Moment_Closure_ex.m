% Analysis of an example gene expression network where the propensity function is not affine (the moment equations are not closed)
% The network incorporates transcriptional feedback from the protein molecules to enhance mRNA degradation. In other words,
% the transcription rate is equal to (k_r + g_fb*[protein number])*[mRNA number].

% Two different moment closure methods are employed, Mass Fluctuation Kinetics(MFK) and Separable Derivative Matching(SDM) and
% are compared with the results obtained from the Stochastic Simulation Algorithm (SSA).
%% 

% Time interval

tspan = [0 10];

% Initial conditions

y0 = [ 0 0 0 0 0];
y0_perturbed = y0 + 10^(-3); %Perturb initial values for SDM integration

%% Mass Fluctuation Kinetics(MFK)  

% Integrating the ODEs

[t_1,Y_1] = ode15s(@first_moment_covariance_MFK_ex, tspan, y0);

%% Separable Derivative Matching(SDM)

% Integrating the ODEs

[t_2,Y_2] = ode15s(@first_moment_covariance_SDM_ex, tspan, y0_perturbed);

%% SSA 

T_final = 10;
P = 10000; % Number of paths
mRNA_SSA=zeros(P,length(t_1));
Protein_SSA=zeros(P,length(t_1));
for k=1:P
[X_results, jumps] = SSA_Moment_Closure_ex(T_final);
    for index=1:length(jumps)
        first_jump_index=find(t_1>jumps(index),1,'first');
        mRNA_SSA(k,first_jump_index:end)=X_results(index+1, 1);
        Protein_SSA(k,first_jump_index:end)=X_results(index+1, 2);
    end
end

%% Moments from SSA

y_SSA=zeros(5,length(t_1));
for index_SSA=1:length(t_1)
    y_SSA(1,index_SSA) = mean(mRNA_SSA(:,index_SSA));
    y_SSA(2,index_SSA) = mean(Protein_SSA(:,index_SSA));
    Sigma_cov = cov(mRNA_SSA(:,index_SSA),Protein_SSA(:,index_SSA));
    y_SSA(3,index_SSA) = Sigma_cov(1,1);
    y_SSA(4,index_SSA) = Sigma_cov(1,2);
    y_SSA(5,index_SSA) = Sigma_cov(2,2);
end
CV_mRNA_SSA = sqrt(y_SSA(3,:))./y_SSA(1,:);
CV_Protein_SSA = sqrt(y_SSA(5,:))./y_SSA(2,:);

%% Plotting MFK - SSA

% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t_1,Y_1(:,1));
hold on
plot(t_1,y_SSA(1,:));
plot(t_1,Y_1(:,2));
plot(t_1,y_SSA(2,:));
hold off
title('First moments of mRNA and protein copy-numbers from MFK')
xlabel('t') 
legend({'y = E_{mRNA} from MFK','y = E_{mRNA} from SSA','y = E_{Protein} from MFK','y = E_{Protein} from SSA'},'Location','northeast')

% Components (distinct) of the covariance matrix
nexttile
plot(t_1,Y_1(:,3));
hold on
plot(t_1,y_SSA(3,:));
plot(t_1,Y_1(:,4));
plot(t_1,y_SSA(4,:));
plot(t_1,Y_1(:,5));
plot(t_1,y_SSA(5,:));
hold off
title('Covariance matrix entries S_{MM}, S_{MP}, S_{PP} from MFK')
xlabel('t') 
legend({'y = S_{MM} from MFK','y = S_{MM} from SSA','y = S_{MP} from MFK','y = S_{MP} from SSA', 'y = S_{PP} from MFK','y = S_{PP} from SSA'},'Location','northeast')

% Coefficient of variation for mRNA and Protein copy-numbers
CV_mRNA_MFK = sqrt(Y_1(:,3))./Y_1(:,1);
CV_Protein_MFK = sqrt(Y_1(:,5))./Y_1(:,2);
nexttile
plot(t_1,CV_mRNA_MFK);
hold on
plot(t_1,CV_mRNA_SSA);
plot(t_1,CV_Protein_MFK);
plot(t_1,CV_Protein_SSA);
hold off
title('Coefficient of variation for mRNA and protein copy-numbers from MFK')
xlabel('t')
legend({'y = CV_{mRNA} from MFK','y = CV_{mRNA} from SSA','y = CV_{Protein} from MFK','y = CV_{Protein} from SSA'},'Location','northeast')

%% Plotting SDM - SSA

% Plotting
tiledlayout(3,1)

% Components of the first moment vector
nexttile
plot(t_2,Y_2(:,1));
hold on
plot(t_1,y_SSA(1,:));
plot(t_2,Y_2(:,2));
plot(t_1,y_SSA(2,:));
hold off
title('First moments of mRNA and protein copy-numbers from SDM')
xlabel('t') 
legend({'y = E_{mRNA} from SDM','y = E_{mRNA} from SSA','y = E_{Protein} from SDM','y = E_{Protein} from SSA'},'Location','northeast')
%
% Components (distinct) of the covariance matrix
nexttile
plot(t_2,Y_2(:,3));
hold on
plot(t_1,y_SSA(3,:));
plot(t_2,Y_2(:,4));
plot(t_1,y_SSA(4,:));
plot(t_2,Y_2(:,5));
plot(t_1,y_SSA(5,:));
hold off
title('Covariance matrix entries S_{MM}, S_{MP}, S_{PP} from SDM')
xlabel('t') 
legend({'y = S_{MM} from SDM','y = S_{MM} from SSA','y = S_{MP} from SDM','y = S_{MP} from SSA', 'y = S_{PP} from SDM','y = S_{PP} from SSA'},'Location','northeast')

% Coefficient of variation for mRNA and Protein copy-numbers
CV_mRNA_SDM = sqrt(Y_2(:,3))./Y_2(:,1);
CV_Protein_SDM = sqrt(Y_2(:,5))./Y_2(:,2);
nexttile
plot(t_2,CV_mRNA_SDM);
hold on
plot(t_1,CV_mRNA_SSA);
plot(t_2,CV_Protein_SDM);
plot(t_1,CV_Protein_SSA);
hold off
title('Coefficient of variation for mRNA and protein copy-numbers from SDM')
xlabel('t')
legend({'y = CV_{mRNA} from SDM','y = CV_{mRNA} from SSA','y = CV_{Protein} from SDM','y = CV_{Protein} from SSA'},'Location','northeast')

%% Estimate coefficients of variation (CV) at stationarity

% CV_mRNA_MFK at stationarity
CV_mRNA_MFK_ss = CV_mRNA_MFK(end,1);

% CV_Protein_MFK at stationarity
CV_Protein_MFK_ss = CV_Protein_MFK(end,1);

% CV_mRNA_SDM at stationarity
CV_mRNA_SDM_ss = CV_mRNA_SDM(end,1);

% CV_Protein_SDM at stationarity
CV_Protein_SDM_ss = CV_Protein_SDM(end,1);

% CV_mRNA_SSA at stationarity
CV_mRNA_SSA_ss = CV_mRNA_SSA(1,end);

% CV_Protein_SSA at stationarity
CV_Protein_SSA_ss = CV_Protein_SSA(1,end);