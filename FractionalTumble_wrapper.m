%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRACKS AND TUMBLES    - Runs track.m                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% TUMBLE_CALC('1K_1522(400).txt');
% function [TUMBLE_MEAN, TUMBLE_STDEV] = TUMBLE_CALC(filename)
posdata = load('1K_1401(400).txt');
maxdisp = 15;
param = struct('dim',2,'quiet',0,'good', 25,'mem', 15);
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
result = track(posdata1, maxdisp, param);
NumParticles = result(end,4);
particleArray = cell(NumParticles,1);
for i = 1:NumParticles
    pdata =  result(result(:,4) == i,:);
    particleArray{i,1} = pdata;
    clear pdata
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TUMBLING RATE                                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for V=1:size(particleArray,1)
%     pdata = particleArray{V,1};
%     VStep = size(pdata,1);
%     for j = 1:VStep-2
%         spx = (pdata(j+2,1) - pdata(j,1))/(pdata(j+2,3) - pdata(j,3));
%         spy = (pdata(j+2,2) - pdata(j,2))/(pdata(j+2,3) - pdata(j,3));
%         vel(j,1) = sqrt(spx^2 + spy^2);
%         if sqrt(spx^2 + spy^2) < 1.388967597
%             TUMBLE(j,1) = 1;
%         else
%             TUMBLE(j,1) = 0;
%         end
%         TUMBLERATE_ALL(V,1) = sum(TUMBLE(:,1) == 1)/length(TUMBLE(:,1));
%     end
% end
% TUMBLE_MEAN = mean(TUMBLERATE_ALL(:,1))
% TUMBLE_STDEV = std(TUMBLERATE_ALL(:,1))/length(TUMBLERATE_ALL(:,1))
% LENGTH = length(TUMBLERATE_ALL(:,1));
counter = 1;
tumble_threshold = 0;
% run_threshold    = 10;

bins = 100; %bins in histogram
edges = linspace(0.1,10,bins+1);
for V=1:size(particleArray,1)
    tumble_num = 1;
    pdata = particleArray{V,1};
    VStep = size(pdata,1);
    for j = 1:VStep-4
        spx_0 = (pdata(j+1,1) - pdata(j,1))/(pdata(j+1,3) - pdata(j,3));
        spy_0 = (pdata(j+1,2) - pdata(j,2))/(pdata(j+1,3) - pdata(j,3));
        spx_1 = (pdata(j+4,1) - pdata(j,1))/(pdata(j+4,3) - pdata(j,3));
        spy_1 = (pdata(j+4,2) - pdata(j,2))/(pdata(j+4,3) - pdata(j,3));
        
        previous_vel = [spx_0 spy_0];
        current_vel = [spx_1 spy_1];
        
        dot_prod = sum(current_vel.*previous_vel);
        mag1 = sqrt(abs(sum(current_vel.*current_vel)));
        mag2 = sqrt(abs(sum(previous_vel.*previous_vel)));
        norm = 1/(mag1*mag2);
        theta = acosd(dot_prod*norm);
    
   if theta > tumble_threshold% & theta_prev < 30
            %mark j as a tumble
            tumble_times(counter,1) = V;
            tumble_times(counter,2) = pdata(j,3);
            %skip j to x frames later
            while pdata(j,3) < tumble_times(counter,2) + 2
                j = j+1;
                if j > length(pdata(j,3))
                    TUMBLE_RATE(V,:) = tumble_num/VStep;
                    break
                end
            end
            tumble_num = tumble_num + 1;
            counter = counter + 1;
        else
            j = j + 1;
   end
  
        theta_prev = theta;
        
    end

   
hold on
    %plot histogram of run times
    if isempty(tumble_times)
        run_lengths = -1;
    elseif length(tumble_times(:,2)) < 2
        run_lengths = -1;
    else
        run_lengths = diff(tumble_times(:,2),1);
    end
    figure(1);
    histogram(run_lengths, edges);
    
end


% [folder, baseFileName, extension] = fileparts(filename);
% statname=[baseFileName '_OUT.txt'];
% save(statname,'TUMBLE_MEAN','TUMBLE_STDEV', 'LENGTH', '-ASCII','-tabs');
% end
