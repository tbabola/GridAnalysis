function ISCstruct = getISCevents(img,fp,name,vidFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    ISCstruct = struct();
    
    %show image and ask for rotation
    figure; imagesc(mean(img,3)); % show z-projection
    x = input('Rotate? (degrees)');
    ISCstruct.rotateDegrees = x;
    if x ~=0
        img = imrotate(img,x);
        imagesc(mean(img,3)); truesize;
    end
    
    [heightImg,widthImg] = size(img,[1 2]);
    sizeSq = 10;

    fprintf('\n Select ISC ROI using the mouse...\n  ')
    roi = drawpolygon;
    ISCstruct.ISCmask = roi; %%save outline to struct
    
    % Grids are organized top to bottom, left to right
    [indices, miniIndices, positiveIndices, miniPosIndices] = getIndices(roi,widthImg,heightImg,sizeSq);
    ISCstruct.indices = indices;
    ISCstruct.miniIndices = miniIndices;
    ISCstruct.posIndices = positiveIndices;
    ISCstruct.miniPosIndices = miniPosIndices;
    
    %plot grids and save
    for i=1:size(positiveIndices,1)
        hold on;
        plot(positiveIndices(i,1:2:end),positiveIndices(i,2:2:end),'Color','g');
    end
    saveas(gcf,[fp '\' name '_gridImage.bmp'])

    %bleachCorrect
    img = bleachCorrect(img,1);
    rois = normalizeGridImg(img,10,positiveIndices);

    [t,n] = size(rois);
    ISCstruct.rois = rois;
    figure; imagesc(corr(rois));
    ISCstruct.corrROIs = corr(rois);
    %ISCstruct.corrshuffROIs = corr(shuffleSignals(rois));

    numToShow = 200;
    if numToShow > size(rois,2)
        numToShow = size(rois,2);
    end
    temp = randperm(n)';
    temp = sort(temp(1:numToShow));
    plotOffset(rois(:,temp),0.2);
    
    %create threshold for ROIs
    roiSignal = ISCstruct.rois;
    medSig = median(roiSignal);
    stdSig = std(roiSignal,1);

    roiThr = roiSignal > medSig + 3*stdSig; %threshold is median + 3 standard deviations

    %
    imgbinary = zeros(ceil(size(img,1)/sizeSq),ceil(size(img,2)/sizeSq),size(rois,1));
    %fill in thresholded ROIs into 3D version to preserve spatial aspects
    %of the signal
    for i = 1:size(miniPosIndices)
        imgbinary(miniIndices(miniPosIndices(i),1),miniIndices(miniPosIndices(i),2),:) = roiThr(:,i);
    end

    imgbinary = imgbinary > 0;
    labels = bwlabeln(imgbinary,26); %looks for connections to any single part of the ROI (vertices, edge, or face)
    count = 0; number = []; wss = [];

    newlabels = zeros(size(labels));
    labelcount = 1;
    for i=1:max(labels,[],'all') 
        [i1,i2,i3] = ind2sub(size(labels),find(labels(:)==i)); %get x,y,t coordinates for label i

        if max(i3)-min(i3) <= 1 | size(i1,1) <= 4 %scrub any events that last less than 3 frames or has less than 4 connected squares
            count = count + 1;
        elseif size(i1,1) > 500
            temp = (labels == i) .* labels;
            singleLabel = labelsToRois(temp, miniIndices,miniPosIndices,size(roiThr));
            generateActivityMovie_labels2(img(:,:,min(i3):max(i3)),singleLabel(min(i3):max(i3),:)',positiveIndices,[],[475 16000]);
            drawnow;

           x = queryEventNum();

           satisfied = 0;
           xscale = 1;
           yscale = 1;
           tscale = 2;
           while ~satisfied
               if x == 1
                   satisfied = 1;
                   newlabels(labels == i) = labelcount;
                   labelcount = labelcount + 1;
               else
                   indices = kmeans([i1/xscale i2/yscale i3/tscale],x); %cluster data into x number of groups
                   temp(sub2ind(size(labels),i1,i2,i3)) = indices; %store these clusters temporarily
                   singleLabel = labelsToRois(temp, miniIndices,miniPosIndices,size(roiThr)); %convert into ROIs for generate activity movie function
                   generateActivityMovie_labels2(img(:,:,min(i3):max(i3)),singleLabel(min(i3):max(i3),:)',positiveIndices,[],[475 16000]);

                   satisfied = querySatisfaction();

                   if ~satisfied
                       [xscale, yscale, tscale] = queryScale(xscale,yscale,tscale);
                   else
                       for j = 1:max(indices)
                            newlabels(temp == j) = labelcount;
                            labelcount = labelcount + 1;
                       end
                   end
               end
           end
        else %if small, just add label into new labels
            newlabels(labels == i) = labelcount;
            labelcount = labelcount + 1;
        end
    end
    disp([num2str(count) ' events were scrubbed.'])
    %%
    labelRoi = labelsToRois(newlabels, miniIndices, miniPosIndices, size(roiThr));
    ISCstruct.labelRoi = labelRoi;
    %%
    if vidFlag
        generateActivityMovie_labels(Kalman_Stack_Filter(single(img)),labelRoi',positiveIndices,[fp '\' name '_ActiveMovie'],[475 max(img,[],'all')*1.5]);
    end
    %% analysis
    ISCstruct.activeArea = sum(sum(labelRoi,1) > 1) / size(labelRoi,2);
    for i = 1:max(newlabels,[],'all')
        [i1,i2,i3] = ind2sub(size(newlabels),find(newlabels(:)==i)); %get x,y,t coordinates for label i
        ISCstruct.event(i).timeStart = min(i3);
        ISCstruct.event(i).timeEnd = max(i3);
        ISCstruct.event(i).eventDuration = max(i3) - min(i3);
        ISCstruct.event(i).area = size(unique([i1 i2],'rows'),1);
        ISCstruct.event(i).maxAmplitude = max(rois(find(labelRoi == i)),[],'all')
    end

    save([fp '\' name '_ISCdata.mat'],'ISCstruct');
end


%helper functions
function [indices, miniIndices, positiveIndices, miniPosIndices] = getIndices(roimask,widthImg,heightImg,sizeSq)
    roiIndices = roimask.Position;
    [indices,miniIndices] = getGrid(widthImg,heightImg,sizeSq);
    positiveIndices = getPositiveGrid(indices,roiIndices);
    [~,miniPosIndices] = ismember(positiveIndices,indices,'rows');
end

function labelRoi = labelsToRois(labels, miniIndices, posIndex, sizeRois)
%LABELROI This function converts label information from 3D grid (x,y,t) back into
%2D list of ROIs (n,t)
    %labels: matrix containing labels for events oriented in 3D space
    %miniIndices: the "mini" grid x,y coordinates
    %posIndex: the "mini" indices for ISC ROIs
    %sizeRois: the size of the final matrix (time,#ISC rois)
    labelRoi = zeros(sizeRois);
    for i = 1:size(posIndex)
        labelRoi(:,i) = squeeze(labels(miniIndices(posIndex(i),1),miniIndices(posIndex(i),2),:));
    end
end

function numEvents = queryEventNum()
%QUERYEVENTNUM This function simply queries how many events the user sees
%with video playback
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
%QUERYSATISFACTION This function simply queries if the user is satisfied
%with the assignment of the ROIs to individual events
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
%QUERYSCALE This function queries the user to change parameters
%(x,y,z-scale) in an attempt to better cluster the events into individual
%clusters. 
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