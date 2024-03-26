function [T, Y] = pend(t_intv,ICs,t_step,m1,m2,L1,L2,gv,method)
    switch method
        case 'rk'
            disp('rk')
            dt = t_intv(1)-t_intv(2)/t_step;
            t = linspace(0,5,t_step);
           
            
            % Preallocate array for solutions
            Y_1 = zeros(4, length(t));
            Y_1(:, 1) = ICs;
            
            % Loop over all time steps
            for i=1:(length(t) - 1)
                k1 = solver(Y_1(:, i), m1, m2, L1, L2, gv);
                k2 = solver(Y_1(:, i) + 0.5*dt*k1, m1, m2, L1, L2, gv);
                k3 = solver(Y_1(:, i) + 0.5*dt*k2, m1, m2, L1, L2, gv);
                k4 = solver(Y_1(:, i) + dt*k3, m1, m2, L1, L2, gv);
                
                Y_1(:, i+1) = Y_1(:, i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
            end
            T = t;
            Y = Y_1.';
        case 'trp'
            disp('trp')
            t0 = t_intv(1);
            t1 = t_intv(2);
            t = linspace(t0,t1,t_step + 1);
           
            h = (t1 - t0) / t_step;
            y = zeros(4, length(t));
            y(:, 1) = ICs;
            
            % Iterate using Trapezoid Method
            for i = 1:t_step
                f = solver(y(:, i), m1, m2, L1, L2, gv);
                Y_pred = y(:, i) + h*f;
                f_next = solver(y(:, i), m1, m2, L1, L2, gv);
                y(:, i+1) = y(:, i) + (h/2)*(f + f_next);
            end
            T = t;
            Y = y.';
        case 'eul'
            disp('eul')
            dt = t_intv(1)-t_intv(2)/t_step;
            t = linspace(0,5,t_step);

            % Preallocate array for solutions
            Y = zeros(4, length(t));
            Y(:, 1) = ICs;

            for i=1:(length(t) - 1)
                % Compute derivatives at current state
                dYdt = solver(Y(:, i), m1, m2, L1, L2, gv);
        
                % Euler's method: y_(n+1) = y_n + h * f(t_n, y_n)
                Y(:, i+1) = Y(:, i) + dt * dYdt;
            end
            Y = Y.';
            T = t;
        otherwise
            disp('Choose rk, trp or eul')
    end
end


function dydt = solver(ICs,m1,m2,L1,L2,g)
    theta1 = ICs(1);
    theta1P = ICs(2);
    theta2 = ICs(3);
    theta2P = ICs(4);
    
    % Compute derivatives
    dydt = zeros(4, 1);
    dydt(1) = theta1P;
    dydt(2) = (-g*(2*m1+m2)*sin(theta1) - m2*g*sin(theta1-2*theta2) - 2*sin(theta1-theta2)*m2*(theta2P^2*L2+theta1P^2*L1*cos(theta1-theta2))) / (L1*(2*m1+m2-m2*cos(2*theta1-2*theta2)));
    dydt(3) = theta2P;
    dydt(4) = (2*sin(theta1-theta2)*(theta1P^2*L1*(m1+m2)+g*(m1+m2)*cos(theta1)+theta2P^2*L2*m2*cos(theta1-theta2))) / (L2*(2*m1+m2-m2*cos(2*theta1-2*theta2)));
end
