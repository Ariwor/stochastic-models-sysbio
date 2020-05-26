% A function that implements the Metropolis-Hastings algorithm to estimate
% the posterior distribution of a parameter 'theta'. For likehood
% computation, the finite state projection (FSP) algorithm is used 
% (function "FSP_MLE.m") for an example auto-repressive gene expression
% network.

function [theta_after] = MH_FSP(theta_prior, sigma)

    % Generate a candidate from the proposal distribution
    
    theta_candidate = lognrnd(log(theta_prior), sigma);
     
    % Compute Hastings ratio
    
    % Assume that the prior distribution is uniform in the range [10, 200]
    Log_likehood_part = FSP_MLE(theta_candidate, 'Data_T.mat') - FSP_MLE(theta_prior, 'Data_T.mat');
    Prior_part = log(unifpdf (theta_candidate,10,200)) - log(unifpdf (theta_prior,10,200));
    Proposal_part = lognrnd(log(theta_candidate), sigma) - lognrnd(log(theta_prior), sigma);
    
    % Calculate the sum
    H =  Log_likehood_part + Prior_part + Proposal_part;
        
    % Convert to linear to calculate the ratio
    H = exp(H);
    
    % Calculate the acceptance probability 
    H = min(1,H);
    
    % Generate a uniform random number (0,1)
    u = rand;
    if u <= H % acceptance
        theta_after = theta_candidate; % set new particle as the candidate
    else % rejection
        theta_after = theta_prior; % set new particle as the previous one
    end
end