%load image
clear all; close all;
addpath(genpath('..\MATLAB Functions'));
[fn, dname] = uigetfile('*.czi'); %prompt for file to analyze
%open image
imgdata = bfopen([dname fn]);
imgdata = imgdata{1,1}; %this is where the imagedata is stored
[fp,name,~] = fileparts([dname fn]); %store filename and path for saving later
img = []; index = 1;
for i = 1:600%1:size(imgdata,1) 
    img(:,:,index) = imgdata{i,1};
    index = index + 1;
end
[fp,name,~] = fileparts([dname fn]);
%%
%ISCstruct = getISCevents(img,fp,name,0);

%%
%IHCstruct = getIHCeventsAuto(img,fp,name,0);

%%
%%make 3D version of signals
bigBin = zeros(size(img),'uint8');
bigISClabel = zeros(size(img),'uint8');
posIndices = ISCstruct.posIndices;
labelROI = ISCstruct.labelRoi > 0;
posRects = [posIndices(:,1) posIndices(:,2) posIndices(:,3)-posIndices(:,1)+1 posIndices(:,6)-posIndices(:,2)+1];
botCircs = [round(IHCstruct.bottomPos) repmat(15,size(IHCstruct.bottomPos,1),1)];

%fill in signals
parfor i=1:size(img,3)
    I = find(labelROI(i,:)==1);
    temp = insertShape(bigBin(:,:,i),'FilledRectangle',posRects(I,:));
    bigBin(:,:,i) = temp(:,:,1);
    
  
    IIHC = find(IHCstruct.thrroi(i,:));
    temp = insertShape(bigBin(:,:,i),'FilledCircle',botCircs(IIHC,:));
    bigBin(:,:,i) = temp(:,:,1);
    bigBin(:,:,i) = bigBin(:,:,i) > 0;
    bigBin(:,:,i) = imgaussfilt(bigBin(:,:,i),[10 10])*255;
end

imgbinary = bigBin > 0;
labels = bwlabeln(imgbinary,26);
count = 0; number = []; wss = [];
centerSqs = [posRects(:,1:2)] + round(posRects(:,3)/2);
newlabels = zeros(size(labels));
labelcount = 1;
temp = zeros(size(labels));

for i=4:max(labels,[],'all') 
    tic;
    [i1,i2,i3] = ind2sub(size(labels),find(labels(:)==i)); %get x,y,t coordinates for label i
    toc;
    if max(i3)-min(i3) <= 1 %scrub any events that last less than 1 frame
        count = count + 1;
    else
        temp = (labels == i); %creates 3D-cluster of just the labeled region;
        singleLabelISCs = labelsToRois(temp, centerSqs, size(ISCstruct.rois)); %3D to 2D ISC rois
        ISClabels = unique((singleLabelISCs > 0).*ISCstruct.labelRoi); %get label numbers for individual ISC events within the 3D blob for seperation of events;
        ISClabels = ISClabels(ISClabels ~=0); %exclude 0s
        centersingleLabelIHCs = labelsToRois(temp, botCircs(:,1:2), size(IHCstruct.thrroi)) > 0; %3D to 2d IHC rois
        
        if size(ISClabels,1) > 1
            satisfied = 0;
            while ~satisfied
               multiLabelISCs = (singleLabelISCs > 0) .* ISCstruct.labelRoi;
               ISClabelsTimes = [];
               for j = 1:size(ISClabels,1)
                   multiLabelISCs(multiLabelISCs == ISClabels(j)) = j;
                   [t1,t2] = find(multiLabelISCs == j);
                   ISClabelsTimes(j,:) = [j min(t1)]; 
               end
               %generateActivityMovie_labels2(img(:,:,min(i3):max(i3)),singleLabel(min(i3):max(i3),:)',posIndices,[],[475 16000]);
               %generateISCIHCActivityMovie(img(:,:,min(i3):max(i3)),multiLabelISCs(min(i3):max(i3),:),posIndices,singleLabelIHCs(min(i3):max(i3),:),IHCstruct.IHCcenters,[],[475 16000])

               x = size(ISClabels,1);

               pkBinary = zeros(max(i3)-min(i3),size(singleLabelIHCs,2));
                    %detect peaks in signal to create points for clustering
                    for j = 1:size(singleLabelIHCs,2)
                         [pks,locs] = findpeaks(IHCstruct.roisignalNorm(min(i3):max(i3),j),'MinPeakProminence',0.3);
                         pkBinary(locs,j) = 1;
                    end
                    [r,c] = find(pkBinary);
                    label = kmeans([r*20 c],x,'MaxIter',200); %time expanded so that this is the major clustering variable, rather than space
                    pkBinaryLab = pkBinary;
                    pkBinaryLab(sub2ind(size(pkBinaryLab),r,c)) = label;

                    figure; imagesc(pkBinaryLab);
                    pause(2);
                    labelTimes = []; mask =[];
                    for j = 1:x
                        [r,c] = find(pkBinaryLab == j);
                        labelTimes(j,:) = [j min(r) max(r)];
                        mask(:,:,j) = pkBinaryLab == j;
                    end
                    labelTimes = sortrows(labelTimes,2);

                    for j= 1:x
                        pkBinaryLab(mask(:,:,labelTimes(j,1))>0) = j;
                    end


                    tempIHCgrouplabels = int16(singleLabelIHCs); %to preserve original structure
                    labelsForEvent = [1];
                    for j = 2:x
                        tempBinary = pkBinaryLab == j;
                        tempBinary = imgaussfilt(single(tempBinary),[0.5 8]) > 0.05;
                        tempBinary = bwlabel(tempBinary);
                        for k=1:max(tempBinary,[],'all')
                            idx = find(tempBinary == k);
                            if k == 1
                                sizeToBeat = size(idx,1);
                                maxIdx = 1;
                            elseif size(idx,1) > sizeToBeat
                                sizeToBeat = size(idx,1);
                                maxIdx = k;
                            end
                        end
                        tempBinary = tempBinary == maxIdx;
                        imagesc(tempBinary);
                        pause(2);
                        timeToChange = ISClabelsTimes(j,2);
                        [r3,c3] = find(tempBinary);
                        groupIHCs = unique([c3]);
                        mask = tempIHCgrouplabels(timeToChange:max(i3),groupIHCs) > 0;
                        tempIHCgrouplabels(timeToChange:max(i3),groupIHCs) = mask * (j);
                        labelsForEvent = [labelsForEvent; j];
                    end
                    generateISCIHCActivityMovie(img(:,:,min(i3):max(i3)),multiLabelISCs(min(i3):max(i3),:),ISCstruct.posIndices,tempIHCgrouplabels(min(i3):max(i3),:),IHCstruct.IHCcenters,[],[0 max(img,[],'all')]); 
                    satisfied = input('Are you satisfied (1) or not (0)?');
                    if ~satisfied %reset original values
                        IHCgrouplabels = tempIHCgrouplabels;
                        oneEvent = input('Is there just one event? (1 or 0)');
                        if oneEvent
                            IHCgrouplabels = singleLabelIHCs;
                            multiLabelISCs = singleLabelISCs;
                            satisfied = 1;
                        else
                            IHCgrouplabels = tempIHCgrouplabels;
                        end
                    else
                        IHCgrouplabels = tempIHCgrouplabels;
                    end
            end   
            
            label3D = RoisToLabels(multiLabelISCs,posRects,IHCgrouplabels,botCircs,size(labels));
            for j = 1:max(label3D,[],'all')
                newlabels(label3D == j) = labelcount;
                labelcount = labelcount + 1;
            end
        else %if small, just add label into new labels
            newlabels(labels == i) = labelcount;
            labelcount = labelcount + 1;
        end
    end
end
disp([num2str(count) ' events were scrubbed.'])

%%
ISClabels = labelsToRois(newlabels, centerSqs, size(ISCstruct.rois));
IHClabels = labelsToRois(newlabels, botCircs(:,1:2), size(IHCstruct.thrroi));

tempISC = ISCstruct.labelRoi > 0;
ISClabels = tempISC .*ISClabels;

tempIHC = IHCstruct.roisignalNorm > IHCstruct.thr;
IHClabels = tempIHC .* IHClabels;

generateISCIHCActivityMovie(Kalman_Stack_Filter(img),ISClabels,posIndices,IHClabels,IHCstruct.IHCcenters,[fp '\' name '_coordNew.mp4'],[475 24000])

%save the coordination
CoorStruct = struct();
CoorStruct.ISClabels = ISClabels;
CoorStruct.posIndices = ISCstruct.posIndices;
CoorStruct.IHClabels = IHClabels;
CoorStruct.IHCcenters = IHCstruct.IHCcenters;
save([fp '\' name '_coorStruct2.mat'],'CoorStruct');

%%

function label3D = RoisToLabels(ISCrois,posRects,IHCrois,botCircs,sizeLabels)
    label3D = zeros(sizeLabels);
    
    maxIHCs = max(IHCrois,[],'all');
    maxISCs = max(ISCrois,[],'all');
    
    for i = 1:maxIHCs
        [r,c] = find(IHCrois == i);
        times = unique(r);
        tempLabel = zeros(sizeLabels);
        for j=1:size(times,1)
            ct = c(r == times(j));
            temp = insertShape(tempLabel(:,:,times(j)),'FilledCircle',botCircs(ct,:));
            tempLabel(:,:,times(j)) = temp(:,:,1);
        end
        label3D(tempLabel>0) = i;
    end
    
    for i = 1:maxISCs
        [r,c] = find(ISCrois == i);
        times = unique(r);
        tempLabel = zeros(sizeLabels);
        for j=1:size(times,1)
            ct = c(r == times(j));
            temp = insertShape(tempLabel(:,:,times(j)),'FilledRectangle',posRects(ct,:));
            tempLabel(:,:,times(j)) = temp(:,:,1);
        end
        label3D(tempLabel>0) = i;
    end
end

function labelRoi = labelsToRois(labels, centerSqs, sizeRois)
    labelRoi = zeros(sizeRois);
    centerSqs(centerSqs > 512) = 512;
    for i = 1:size(centerSqs,1)
        labelRoi(:,i) = squeeze(labels(centerSqs(i,2),centerSqs(i,1),:));
    end
end

function numEvents = queryEventNum()
    x = input('How many events do you see? ');
    
       while 1
           if isnumeric(x) && x > 0
               numEvents = x;
               return
           else
             x = input('How many events do you see? (Enter a postive number) ');
           end
       end
end

function satisfied = querySatisfaction()
    satisfied = input('Satisfied?');
    
    while 1
       if isnumeric(satisfied) && satisfied >= 0 && satisfied <=1
           return
       else
         satisfied = input('Are you satisfied? (0) No (1) Yes');
       end
    end
end

function [xscale,yscale,tscale] = queryScale(xscale,yscale,tscale)
    xscale = xscale; yscale = yscale; tscale = tscale;
    x = input(['Which scale would you like to alter to adjust clustering?\n' ... 
        '(1) x, Current value: ' num2str(xscale) '\n' ...
        '(2) y, Current value: ' num2str(yscale) '\n' ...
        '(3) t, Current value: ' num2str(tscale) '\n']);
    
    while 1
       if x > 0 && x <= 3
           if x == 1
                xscale = input('What is the new scale? ');
                while 1
                     if isnumeric(xscale)
                           return
                     else
                        xscale = input('What is the new scale? (Must be numeric) ');
                     end
                end
           elseif x == 2
                yscale = input('What is the new scale?');
                while 1
                     if isnumeric(yscale)
                           return
                     else
                        yscale = input('What is the new scale? (Must be numeric) ');
                     end
                end
           elseif x == 3
                tscale = input('What is the new scale?');
                while 1
                     if isnumeric(tscale)
                           return
                     else
                        tscale = input('What is the new scale? (Must be numeric) ');
                     end
                end
           end
       else
         x = input(['Which scale would you like to alter to adjust clustering?\n' ... 
        '(1) x, Current value: ' num2str(xscale) '\n' ...
        '(2) y, Current value: ' num2str(yscale) '\n' ...
        '(3) t, Current value: ' num2str(tscale) '\n']);
       end
   end
end

%%
%%unused clustering code
%                x = queryEventNum();
%                
%                satisfied = 0;
%                xscale = 1/10;
%                yscale = 1/10;
%                tscale = .01;
%                while ~satisfied
%                    if x == 1
%                        satisfied = 1;
%                        newlabels(labels == i) = labelcount;
%                        labelcount = labelcount + 1;
%                    else
%                        indices = kmeans([i1/xscale i2/yscale i3/tscale],x); %cluster data into x number of groups
%                        temp(sub2ind(size(labels),i1,i2,i3)) = indices; %store these clusters temporarily
%                        
%                        startTimes =[];
%                        for j=1:size(ISClabels,1)
%                            [t1,t2] = find(multiLabelISCs == j);
%                            if ~isempty(t1)
%                                startTimes(j,:) = [j min(t1)];
%                            else
%                                startTimes(j,:) = [j 0];
%                            end
%                        end
%                        
%                        [~,sortI] = sort(startTimes(:,2));
%                        startTimes = startTimes(sortI,:);
%                        
%                        for j=2:x
%                            indices(i3 <= startTimes(j,2) & indices == startTimes(j,1)) = startTimes((j-1),1);
%                        end
%                        temp(sub2ind(size(labels),i1,i2,i3)) = indices; %store these clusters temporarily
%                        tempIHC = labelsToRois(temp, botCircs(:,1:2), size(IHCstruct.thrroi));
%                        
%                        generateISCIHCActivityMovie(img(:,:,min(i3):max(i3)),multiLabelISCs(min(i3):max(i3),:),posIndices,tempIHC(min(i3):max(i3),:),IHCstruct.IHCcenters,[],[475 16000])
% 
%                        %satisfied = 1;
%                        satisfied = querySatisfaction();
% 
%                        if ~satisfied
%                            [xscale, yscale, tscale] = queryScale(xscale,yscale,tscale);
%                        else
%                            for j = 1:max(indices)
%                                 newlabels(temp == j) = labelcount;
%                                 labelcount = labelcount + 1;
%                            end
%                        end
%                    end
%                end