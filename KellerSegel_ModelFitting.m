% Keller-Segel Model in a Linear Chemoattractant Concentration
% B_t = {mu * B_x - V * B}_x
% V = 8v/3pi * tanh {[chi_0*pi/8v])* [k_d(k_D+C)] * C_x}
% Bacterial Concentration: B(x,t)
% Chemoattractant Concentration: C(x)=C_0*x/W
% Simulation for a 1mm Channel across a 15 minute simulation

clear; close all
tic
x = 0:1E-6:1E-3;
t = 0:60:3600;
%model determined by pdede solver
P0 = [10E-6, 10E-8, 5.9E-10];
sol = kellersegel(x,t,P0);
B = sol(:,:,1);
C= 100*exp(-x);
%Adding gaussian noise to signal
sol_noise = awgn(sol,10,'measured');
%Modified output with noise i.e. data to fit
B_n=B;

%initial guess at parameters
% p1 = [20E-6, 10E-8,5E-10]
%function to optimise
% fun = @(p1)sseval(x,t,p1);
% options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');

% P = lsqnonlin(fun,p1, [],[], options)

%plotting model against data
% model = kellersegel(x,t,P);
% B_m = model(:,:,1);
% plot(x,B_n(900,:),'r-', x, B_m(900,:),'b-')
% legend("modified signal", "model");

%VARIABLES TO FIT
% VEL - P(1)
% CHI - P(2)
% MU  - P(3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTTING OF B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
subplot(1,2,1);
% dist=linspace(0,1,1000);
% conc=exp(-1000.*abs(dist));
plot(x,B(1,:),'Color','r','LineW',1) % Intial Condition
hold on
plot(x,B([60],:),'k', 'LineW',1.5)
plot(x,C, 'Color', [.8 .8 .8], 'LineW',1)
% plot(x,Bn([900],:),'b', 'LineW', 0.8)
ylim([0 2]);
legend('B(x,t=0): Intial Condition', 'B(x,t=900)','Location', 'northwest')
xlabel('Displacement, x, [m]'); ylabel('Concentration, [mM]')
hold off
subplot(1,2,2);
s = surf(x,t,sol);
colormap winter
c = colorbar('Location', 'Northoutside');
c.Label.String = 'Concentration [mM]';
set(s,'LineStyle','none')
xlabel('Displacement, x, [m]'); ylabel('Time, t, [s]');
zlabel('Concentration of Bacteria, B(x,t), [mM]');
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sse = sseval(x,t,P)
%sum of square of residuals
P0 = [30E-6, 10E-8, 5.9E-10];
sol = kellersegel(x,t,P0);
sol_noise = sol;% awgn(sol,10,'measured');
 sse = sum((sol_noise - kellersegel(x,t,P)).^2);
end

function [sol]= kellersegel(x,t,P0)
% ROUTINE WHICH NUMERICALlY SOLVES KELLER-SEGEL EQN USING PDEPE
sol = pdepe(0,@pdecomp,@intial,@boundary,x,t);
    function [c,f,s] = pdecomp(x,t,B,DBDx)
    % PDE COMPONENTS
    C0 = 100;           % Chemoattractant Initial Concentration: mM
    W = 12E-2;           % Width of Chanell: metres, m
%     C = C0*x/W;         % Chemoattractant Concentration: mM
C= C0*exp(-x);
%     P0(3) = 5.9E-10;      % Random Motility Coefficient: m^2 / s
%     P0(2) = 5E-8;      % Chemotatic Coefficient: m^2 / s
%     P0(1) = 30E-6;     % Average Velocity of Bacteria m/s
    kd=0.125;             % Dissociation Constant: mM
 
    Vd = (8*P0(1)/(3*pi))*tanh((P0(2)*pi*kd/(8*P0(1)*(kd+C)^2))*C0/W);

    c = 1;
    f = P0(3)*DBDx - Vd*B;
    s = 0;
    end

    function [pl,ql,pr,qr] = boundary(xl,Bl,xr,Br,t)
    % BOUNDARY CONDITIONS
    pl = 0; pr = 0; 
    ql = 1; qr = 1;
    end

    function B0 = intial(x)
    % INTIAL CONDITIONS
    B0 = 1;
    end
end