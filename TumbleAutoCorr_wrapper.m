

close all
% pnum = 28; % trajecotry number
for pnum=1:182
clear d;
clear phi;
clear phi_b;
clear a;
a = particleArray{pnum,1}; % load the trajectory
% Nkeep = 100;
% a=b(1:Nkeep,:);
% disp(length(a)); % dispay to the command window what the length of the trajecotry is in frames.

% plot the trajectory
% plot(a(:,1),a(:,2), 'k');
% hold on;
% % tester to see if detected locations do indeed have a tumble.
% scatter(a(4:8,1),a(4:8,2), 'r');
% scatter(a(41:53,1),a(41:53,2), 'r');
% scatter(a(58:60,1),a(58:60,2), 'r');
% scatter(a(63:71,1),a(63:71,2), 'r');
% scatter(a(79:82,1),a(79:82,2), 'r');

%%

% clear the figure
% hold off;

% get the displacement vectors of each of the positions along the
% trajectory.
d(:,1) = a(2:end,1) - a(1:end-1,1);
d(:,2) = a(2:end,2) - a(1:end-1,2);
d(:,3) = sqrt( (a(2:end,1) - a(1:end-1,1)).^2 + (a(2:end,2) - a(1:end-1,2)).^2 );

N = 4; % number of time steps over which you expect the tumble to last (smaller better for short trjectories))
cons = zeros(length(d)-1,1);
phi = zeros(length(d)-N-1,1);

% loop to get the consecutive dot products of the direction vectors
% normalised by their lengths.
for i = 1:length(d) - 1
    
        cons(i) = abs( (d(i,1) * d(i+1,1)) + (d(i,2) * d(i+1,2)) )/ (sqrt(d(i,1)^2 + d(i,2)^2)*sqrt(d(i+1,1)^2 + d(i+1,2)^2));
    
end

% get the order parameter, check this: it doesn't normalize correctly
% currently, unclear why.
for i = 1:length(cons) - N
    
    phi(i,:) = 1/(N+1) * sum(cons(i:i+N));
    
end

% determine a threshold for the order parameter to decide where motion is
% no longer balistic.
thresh = 0.5;

% binarise the trajectory (1 = balistic, 0 = non-(.))
for i = 1:length(phi)
    
    if(phi(i) >= thresh)
        phi_b(i,:) = 1;
    else
        phi_b(i,:) = 0;
    end
    
end
average_v = mean(d(1:end-5,3))* 0.1799898*25;
average_phi = mean(phi);
% plot the binarised order parameter to show where tumbles have been
% detected.
plot(d(1:end-5,3)* 0.1799898*25, phi, '.','Color', 'k')
% plot (average_v, phi_b, '.', 'Color', 'k')
% plot(phi, 'LineWidth', 1);
hold on
% plot(d(:,3)* 0.1799898*25);
box on;
hXLabel = xlabel('Velocity');
hYLabel = ylabel('{\psi(N)}');
set([hXLabel, hYLabel], 'FontName', 'CMU Serif')
set([hXLabel, hYLabel], 'FontSize', 24);
set(gca, 'FontName', 'CMU Serif');
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',20)
set(gca,'XTickLabelMode','auto')
b= get(gca,'YTickLabel');  
set(gca,'YTickLabel',b,'fontsize',20)
set(gca,'YTickLabelMode','auto')
% xlim([0 40])
ylim([-0.1 1.1]);
set(gca, 'Box', 'on', 'TickDir', 'in', 'TickLength', [.01 .01], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'off', ...
    'XColor', [0 0 0], 'YColor', [0 0 0], 'LineWidth', 1);
end