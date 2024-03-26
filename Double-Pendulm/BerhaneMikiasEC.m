 close all;
%% Double Pendulum System Parameters
gv = 9.81; % Gravity Constant
m1 = 1;    % Mass of bob 1 
m2 = 1;    % Mass of bob 2
L1 = 1;    % Length of Rod 1
L2 = 1;  % Length of Rod 2 (Cases: [1/4 1/2] m)

%% Initial Conditions
ICs = [pi/3 0 pi*2/3 0];
y = [pi/3 0 pi*2/3 0];
t = 500;

%% Time Step 
t_intv = [0 5];
t_step = 500;

%% Integrators

[t_rk,y_rk]   = pend(t_intv,ICs,t_step,m1,m2,L1,L2,gv,'rk');
[t_trp,y_trp] = pend(t_intv,ICs,t_step,m1,m2,L1,L2,gv,'trp');
[t_eul,y_eul] = pend(t_intv,ICs,t_step,m1,m2,L1,L2,gv,'eul');

%% ODE 45 

[T,Y] = ode45(@(t,y) YDOT(t,y,gv,m1,m2,L1,L2),t_intv,ICs);

%% Movie
movie(Y, L1, L2, "Ode45")
movie(y_rk, L1, L2, "Runge Kutta")
movie(y_trp, L1, L2, "Trapezoid")
movie(y_eul, L1, L2, "Euler")

%% Plot
hold on
plot(t_eul,y_eul(:,1),'--','color',"#EDB120",'LineWidth',2)
plot(t_rk,y_rk(:,1),':','color',"#0072BD",'LineWidth',2)
plot(t_trp,y_trp(:,1),'-.','color',"#D95319",'LineWidth',2)
plot(T,Y(:,1),'k-','LineWidth',2)
xlabel ('Time [s]')
ylabel ('\theta_1 [rad]')
legend('Euler','Trapezoid','Runge-Kutta','ODE45','Location','best')
ylim([-5 5])
grid on
box on
hold off