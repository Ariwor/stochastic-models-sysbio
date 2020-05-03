% A system of ODEs that includes 2 equations for the first moment vector and 3 for the
% 3 distinct entries of the covariance matrix of a simple gene expression network,
% that also incorporates transcriptional feedback from the protein molecules. In other
% words, the transcription rate is equal to k_0 + k+1*[protein number].

% In this case though, k_1=+1 and as a result there is positive feedback

function dydt = first_moment_covariance_matrix_pos(t,y)
k_0=10;
k_1=1;
k_p=5;
g_r=10;
g_p=2;
 dydt = zeros(5,1);
 dydt(1) = -g_r*y(1) + k_1*y(2) + k_0; % dM/dt
 dydt(2) = k_p*y(1) - g_p*y(2); % dP/dt
 dydt(3) = -2*g_r*y(3) + 2*k_1*y(4) + k_1*y(2) + g_r*y(1) + k_0; % dS_MM/dt
 dydt(4) = k_1*y(5) + k_p*y(3) - (g_r + g_p)*y(4); % dS_MP/dt
 dydt(5) = -2*g_p*y(5) + 2*k_p*y(4) + k_p*y(1) + g_p*y(2); % dS_PP/dt
end

