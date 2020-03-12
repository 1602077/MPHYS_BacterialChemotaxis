%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRACKS AND PLOTS DATA - Runs track.m                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
% posdata = load('1K_1450(400).txt'); #28 contains track used in report
posdata = load('4K_1044(400).txt');
maxdisp = 15; 
param = struct('dim',2,'quiet',0,'good', 25,'mem', 15);
% dim:   Dimensions of System
% quiet: Removes output text
% good:  # of frames particles need to be present for, else thrown away
% mem:   memory, # of frames a particle can dissapear for.
% pos:   array of form [x,y,t]

keep = ones(size(posdata,1),1);

for i = 1:size(keep,1)
    if posdata(i,1) ~= 0
        continue
    else
        if posdata(i,2) == 0 && posdata(i,3) == 0
            keep(i,1) = 0;
        end
    end
end

posdata1 = posdata(keep == 1,:);
% result: array of form [x,y,t,id] - id = particle id #
% result = track(posdata1,maxdisp,param);
result = track(posdata1, maxdisp, param);
NumParticles = result(end,4);
particleArray = cell(NumParticles,1);

for i = 1:NumParticles
    % Select all the rows from result such that the
    % particle number is the current index value
    pdata =  result(result(:,4) == i,:);
    particleArray{i,1} = pdata;
    %i is the particle id number
    clear pdata    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTTTING TRAJECTORIES                                             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fh = figure;
% set(fh,'color','white'); box on; hold on;
% xlabel('x displacement, pix'); ylabel('y displacement, pix');
% % title('25/11/2019 spec. 00ug/mol 11:52AM (1000 Frames');
% hold on;
% CM = jet(2*NumParticles);
% maxX = 2048;
% maxY = 2048;
% xlim([0 2048]);
% ylim([0 2048]);
% 
% for i =1:NumParticles
% 
%     pdata = particleArray{i,1};
% %     plot(pdata(:,1),pdata(:,2),'color','k');
%     hold on;
%     text(pdata(1,1),pdata(1,2),num2str(i)); hold on;
%     maxX = max([maxX;pdata(:,1)]);
%     maxY = max([maxY;pdata(:,2)]);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATING AVERAGE VELOCITY                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for V=1:size(particleArray,1) %choosing correleation between two specific particles
pdata = particleArray{V,1};
VStep = size(pdata,1);
    for j = 1:VStep-1
        spx = (pdata(j+1,1) - pdata(j,1))/(pdata(j+1,3) - pdata(j,3));
        spy = (pdata(j+1,2) - pdata(j,2))/(pdata(j+1,3) - pdata(j,3));
        vel(j,1) = sqrt(spx^2 + spy^2);
        myvel = vel; 
%  myvel called by autocorrelation, required as clear vel on each iteration
    end
sdvel=std(vel(:))/ sqrt(length(vel));
newvel=vel(vel(:)>(mean(vel(:))-sdvel));
Avel(V,:)=mean(newvel);
 clear vel  
end
%previsouly Avel
Average=mean(Avel);
VELN = Avel * 0.1799898*25
Vel_Out = Average * 0.1799898*25 % micrometres /s
sd_out = sdvel * 0.1799898*25

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TUMBLING RATE                                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for V=1:size(particleArray,1)
% pdata = particleArray{V,1};
% VStep = size(pdata,1);
%     for j = 1:VStep-5
%         spx = (pdata(j+5,1) - pdata(j,1))/(pdata(j+5,3) - pdata(j,3));
%         spy = (pdata(j+5,2) - pdata(j,2))/(pdata(j+5,3) - pdata(j,3));
%         vel(j,1) = sqrt(spx^2 + spy^2);
%         if sqrt(spx^2 + spy^2) < 1.388967597
%             TUMBLE(j,1) = 1;
%         else
%             TUMBLE(j,1) = 0;
% %         RATE = sum(TUMBLE(:,1) == 1)/length(TUMBLE(:,1))
%         end
% %         RATE = 
% TUMBLERATE_ALL(V,1) = sum(TUMBLE(:,1) == 1)/length(TUMBLE(:,1));
%     end
% end
% TUMBLE_MEAN = mean(TUMBLERATE_ALL(:,1))
