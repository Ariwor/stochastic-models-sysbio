%% Moment Dynamics
function dydt = first_second_moments(t, y)
 kappa = 2;
 gamma = 1;
 dydt = zeros(2,1);
 dydt(1) = kappa - gamma*y(1);
 dydt(2) = -2*gamma*y(2) + (2*kappa + gamma)*y(1) + kappa;
end
