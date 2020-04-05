function IHCstruct = getIHCeventsAuto(img,fp,name,vidFlag)
%getIHCeventsAuto Automatically detects rois and generates a structure
%containing useful analysis of signals
%   Detailed explanation goes here

    meanImg = mean(img,3);
    meanImg = (meanImg - mean(meanImg,'all'))/std(meanImg,[],'all'); %image is normalized

    figure(1); hold off; imagesc(meanImg);
    fprintf('\n Select IHC search area using the mouse...\n  ')
    roi = drawpolygon;
    IHCmask = createMask(roi);
    IHCstruct = struct();
    IHCstruct.IHCmask = IHCmask;
    IHCstruct.meanImg = meanImg;

    %enhance nuclei darkness before continuing
    enhSatisfied = 0;
    while ~enhSatisfied
        meanImgEnh = enhanceNuclei(meanImg);
        maskedImg = meanImg.*IHCmask;
        imagesc(imbinarize(imcomplement(maskedImg),'adaptive','Sensitivity',1));
        enhSatisfied = input('Are you satisfied (1) or would you like to reset (0)?');
    end
    IHCstruct.meanImgEnh = meanImgEnh;
    meanImg = meanImgEnh;
    
    %
    %thresholding method
    figure(1); 
    
    binarized = imbinarize(imcomplement(maskedImg),'adaptive','Sensitivity',1);
    structs = regionprops(binarized,'Centroid','Area','PixelList');
    structs = structs([structs.Area] < 200 & [structs.Area] > 5); %remove very small and very large detected regions
    
    roiHandles = {};
    figure(2); imagesc(meanImg);
    for i=1:size(structs,1)
        roiHandles{i} = drawpoint('Position',structs(i).Centroid);
    end
    
    %add any points manually
    disp('Add ROIs manually, then hit escape.')
    points = []; roi = []; addedPoints = 0;
    while 1
        temp = drawpoint;
        if isempty(temp.Position)
            break
        else
            roiHandles{end+1} = temp;
            addedPoints = 1;
        end
    end

    for i=1:size(roiHandles,2)
        points(i,:) = roiHandles{i}.Position;
    end

    if addedPoints
        tempStr = input('Sort by x or y?','s');
        if strcmp(tempStr,'x')
            [~,temp] = sort(points(:,1));
            points = points(temp,:);
        elseif strcmp(tempStr,'y')
            [~,temp] = sort(points(:,2));
            points = points(temp,:);
        end
    end
    
    
    %% fit a curve to the points, then drag line to medial side of IHCs to move centers towards the medial edges for ROI drawing
    f = fit(points(:,1),points(:,2),'poly4');
    tempStr = [min(points(:,1)):0.5:max(points(:,1))]';
    polypoints = [tempStr f(tempStr)];

    roi = drawpolyline('Position',polypoints(16:15:end-15,:));
    input('Drag line to the medial side of the IHCs, then press enter.')
    newpolypoints = polypoints - (polypoints(16,:) - roi.Position(1,:));
    IHCstruct.polypoints = newpolypoints;
    
    figure(2); imagesc(meanImg);
    bottomPos = [];
    distTowardsLine = 7; %number of pixels to move towards the line, adjust this to move location of the roi

    for i = 1:size(points,1)
          %go from centroid towards drawn line
          tempCentroid = points(i,:);
          dist = pdist2(tempCentroid,newpolypoints);
          [~,I]= min(dist);
          tempPoint = newpolypoints(I,:);
          tempslope = tempCentroid - tempPoint;
          %a little geometry
          tempSq = tempslope.^2;
          tempSqAdd = sum(tempSq);
          unitToMove = sqrt(distTowardsLine^2/tempSqAdd);
          tempPoint = tempCentroid - tempslope.*unitToMove;
          bottomPos(i,:) = tempPoint;
    end
    
    %uncomment for troubleshooting
%     figure(2); imagesc(meanImg);
%     for i = 1:size(bottomPos,1)
%         drawpoint('Position',points(i,:),'Color','r');
%         drawpoint('Position',bottomPos(i,:),'Color','g');
%     end
    
    %plot ellipses and confirm with user that these are appropriately
    %placed
    
    I = {}; roihandles = {};
    figure(2); imagesc(meanImg); truesize; 
    for i=1:size(bottomPos,1)
        roihandles{i} = drawellipse('Center',bottomPos(i,:),'SemiAxes',[3 3]);
    end

    input('Adjust any points and then hit enter at the command line.');
    
    count = 1;
    for i=1:size(points,1)
        h = roihandles{i};
        if ishandle(h)
            finalPoints(count,:) = h.Position;
            mask = h.createMask;
            I{count} = find(mask);
            count = count + 1;
        else
            points(i,:) = [];
        end
    end
    IHCstruct.IHCcenters = points;
    
    roisignal = extractROIs(img, I);
    
    %normalize and smooth rois
    prcF = prctile(roisignal,5); %5th percentile as Fo
    roisignalNorm = (roisignal - prcF) ./prcF;
    for i = 1:size(roisignalNorm,2)
        roisignalNorm(:,i) = smooth(roisignalNorm(:,i),3);
    end
    
    medians = []; stds = []; 
    window = 15; %number of frames to analyze
    %take windows of measurements for median and std for more accurate
    %measurements (i.e. prevents active cells from being punished)
    for i = 1:size(roisignalNorm,1)-window
       medians(i,:) =  median(roisignalNorm(i:i+window-1,:));
       stds(i,:) = std(roisignalNorm(i:i+window-1,:),[],1);
    end

    finMedian = median(medians);
    finStds = median(stds);
    thr = finMedian + 4*finStds;

    binary = roisignalNorm > thr;
    
    %generate movie 
    if vidFlag
        generateActivityMovieIHCs(Kalman_Stack_Filter(single(img)), binary, points,[fp '\' name '_IHCactivityMovie'],[min(img,[],'all') 20000]);
    end
    
    IHCstruct.bottomPos = finalPoints;
    IHCstruct.roisignal = roisignal;
    IHCstruct.roisignalNorm = roisignalNorm;
    IHCstruct.roiMedians = finMedian;
    IHCstruct.roiStds = finStds;
    IHCstruct.thr = thr;
    IHCstruct.thrroi = binary;
    save([fp '\' name '_IHCstruct.mat'],'IHCstruct');

end

function meanImgEnh = enhanceNuclei(meanImg)
%MEANIMGENH This function increases the contrast of the nuclei. Prompts to click on nucleus, and then decreases the intensity of pixel values below the selected pixel value
%   meanImg: mean image across the entire stack
    figure(1);
    meanImgEnh = meanImg;
    disp('Click within a nuclei to enhance contrast. Press escape to end enhancement.')
    while 1
        hold off; imagesc(meanImgEnh);
        roi = drawpoint();
        if isempty(roi.Position)
            break
        end
        pos = roi.Position; 
        val = meanImg(round(pos(1)),round(pos(2)))
        tempImg = meanImgEnh;
        tempImg(tempImg < val) = tempImg(tempImg < val) * 0.9;
        imagesc(tempImg);
        meanImgEnh = tempImg;
    end
end

function angles = getAngles(bottomPos)
%GETANGLES Returns the angle of the hair cell on the curve. This would be
%useful if you wanted to rotate the ROI by a certain amount. Not used in
%this particular version but it is an option.

     movAvg = 5;
     slope = 1;
     points = bottomPos;
     points(:,2) = 512 - points(:,2) ;
     
     for i=1:size(points,1)
        startPt = i - floor(movAvg/2);
        endPt = i + floor(movAvg/2);
        startPt;
        if startPt < 1
            startPt = 1;
        elseif endPt > size(points,1)
            endPt = size(points,1);
        end
        p = polyfit(points(startPt:endPt,1),points(startPt:endPt,2),1);
        slope = [slope; p(1)];
    end
    
    angles = atan(slope) * 180/pi;
    angles = angles + 90;
end

function roisignal = extractROIs(img, roiMasks)
    roisignal = zeros(size(img,3),size(roiMasks,1)); %will hold masked ROIs across the image
    %go through image and extract rois
    for i=1:size(img,3)
        temp = img(:,:,i);
        for j=1:size(roiMasks,2)
            tempMask = roiMasks{j};
            roisignal(i,j) = mean(temp(tempMask),'all');
        end 
    end
end