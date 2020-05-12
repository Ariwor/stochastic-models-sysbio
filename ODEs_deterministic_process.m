% A system of ODEs for the deterministic concentration process x(t) = (mRNA(t), Protein(t)) for an example gene expression network,
% which incorporates transcriptional feedback from the protein molecules to enhance mRNA degradation.

function dxdt = ODEs_deterministic_process(t,x)

 k_r = 1;
 g_r = 10;
 k_p = 5;
 g_p = 2;
 g_fb = 6;

 dxdt = zeros(2,1);
 dxdt(1) = k_r -(g_r + g_fb*x(2))*x(1); % dM/dt
 dxdt(2) = k_p*x(1) - g_p*x(2); % dP/dt
end

