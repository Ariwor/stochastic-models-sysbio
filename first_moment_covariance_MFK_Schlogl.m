% A system of ODEs that includes 2 equations for the first moment vector and 3 for the
% 3 distinct entries of the covariance matrix for the Schlogl model.
% Rate constants depend on a volume parameter V, which is defined by the user.

% Mass Fluctuation Kinetics(MFK) is used to "close" the moment equations.

function dydt = first_moment_covariance_MFK_Schlogl(t,y,V)

 c_1 = V;
 c_2 = 1;
 c_3 = 5/V;
 c_4 = 0.2/V;
 c_5 = 5;

 dydt = zeros(5,1);
 
 S_112 = 0;
 S_221 = 0;
 S_111 = 0;
 
 dydt(1) = c_1 + 2*c_2*y(2) - (c_3/2)*(y(1)^(2) - y(1) + y(3)) - c_4*(y(1)*y(2) + y(4)) - c_5*y(1); % dE(X)/dt
 dydt(2) = -c_2*y(2) + (c_3/2)*(y(1)^(2) - y(1) + y(3)); % dE(Y)/dt
 dydt(3) = 2*y(3)*(-c_3*y(1) + c_3/2 - c_4*y(2) - c_5) + 2*y(4)*(2*c_2 - c_4*y(1)) + c_1 + 4*c_2*y(2) + (c_3/2)*(y(1)^(2) - y(1) + y(3)) + c_4*(y(1)*y(2) + y(4)) + c_5*y(1) - c_3*S_111 - 2*c_4*S_112; % dS_XX/dt
 dydt(4) = y(4)*(-c_3*y(1) + c_3/2 - c_4*y(2) - c_5) + y(5)*(2*c_2 - c_4*y(1)) + y(3)*(c_3*y(1) - c_3/2) - y(4)*c_2 - 2*c_2*y(2) - (c_3/2)*(y(1)^(2) - y(1) + y(3)) + (c_3/2)*(S_111 + S_112) - c_4*S_221; % dS_XY/dt
 dydt(5) = 2*y(4)*(c_3*y(1) - c_3/2) - 2*y(5)*c_2 + c_2*y(2) + (c_3/2)*(y(1)^(2) - y(1) + y(3)) + c_3*S_112; % dS_YY/dt
end

