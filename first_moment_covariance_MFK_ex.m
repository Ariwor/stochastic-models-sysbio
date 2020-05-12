% A system of ODEs that includes 2 equations for the first moment vector and 3 for the
% 3 distinct entries of the covariance matrix of an example gene expression network,
% which incorporates transcriptional feedback from the protein molecules to enhance mRNA degradation.
% In other words, the transcription rate is equal to (k_r + g_fb*[protein number])*[mRNA number].

% Mass Fluctuation Kinetics(MFK) is used to "close" the moment equations.

function dydt = first_moment_covariance_MFK_ex(t,y)
 k_r=1;
 g_r=10;
 k_p=5;
 g_p=2;
 g_fb=6;
 
 dydt = zeros(5,1);
 
 S_112 = 0;
 S_221 = 0;
 
 dydt(1) = k_r - g_fb*y(4) - (g_r + g_fb*y(2))*y(1); % dE(mRNA)/dt
 dydt(2) = k_p*y(1) - g_p*y(2); % dE(Protein)/dt
 dydt(3) = -2*y(3)*(g_r + g_fb*y(2)) - 2*y(4)*g_fb*y(1) + k_r + (g_r + g_fb*y(2))*y(1) + g_fb*y(4) - 2*g_fb*S_112; % dS_MM/dt
 dydt(4) = -y(4)*(g_r + g_fb*y(2) + g_p) + y(3)*k_p - y(5)*g_fb*y(1) - S_221*g_fb; % dS_MP/dt
 dydt(5) = 2*y(4)*k_p - 2*y(5)*g_p + k_p*y(1) - g_p*y(2); % dS_PP/dt
end


