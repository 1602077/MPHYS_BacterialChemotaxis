%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRETRACKING ROUNTINES - Runs bpass.m, pkfnd.m, cntrd.m for i frames %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear
% i: Number of frames recorded
delete 1_1016(1-400).txt
for i = 1:400
    if i<10
        FileName = sprintf('1_1016_T00%d.tif', i);
    elseif i>=10 & i<100
        FileName = sprintf('1_1016_T0%d.tif', i);
    elseif i>=100 & i>400
        FileName = sprintf('1_1016_T%d.tif', i);
    end  
    if exist(FileName, 'file')
%         figure();
        a = double(imread(FileName));
%         colormap('gray'),imagesc(a);
        aprime = imcomplement(a);
        b=bpass(aprime,2,10);
%         imshow(b)
        pk = pkfnd(b,1E3,5);
        cnt = cntrd(b,pk,15);
%         whos cnt;
%         max(max(b))
%         hist(mod(cnt(:,1),1),10);
        x = cnt(:,1);
        y = cnt(:,2);
        L = length(x);
        % Printing frame number to corresponding [x y] value of centroid
        I = repmat(i, [1, L]); 
        fileID = fopen('1_1016(1-400).txt', 'a');
        fprintf(fileID, '%f\t%f\t%f\n', [x,y, I.']');
        fclose(fileID);
    else
        fprintf('File %s does not exist.\n', FileName);
    end
end


