function [dydt] = YDOT(t, y, g, m1, m2, L1, L2)

    theta1 = y(1);
    theta1P = y(2);
    theta2 = y(3);
    theta2P = y(4);
    
    % Compute derivatives
    dydt = zeros(4, 1);
    dydt(1) = theta1P;
    dydt(2) = (-g*(2*m1+m2)*sin(theta1) - m2*g*sin(theta1-2*theta2) - 2*sin(theta1-theta2)*m2*(theta2P^2*L2+theta1P^2*L1*cos(theta1-theta2))) / (L1*(2*m1+m2-m2*cos(2*theta1-2*theta2)));
    dydt(3) = theta2P;
    dydt(4) = (2*sin(theta1-theta2)*(theta1P^2*L1*(m1+m2)+g*(m1+m2)*cos(theta1)+theta2P^2*L2*m2*cos(theta1-theta2))) / (L2*(2*m1+m2-m2*cos(2*theta1-2*theta2)));
end
