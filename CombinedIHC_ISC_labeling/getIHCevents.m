function IHCstruct = getIHCeventsManual(img,fp,name)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    figure; imagesc(std(img,[],3));

    points = []; roi = [];
    while 1
        temp = drawpoint;
        if isempty(temp.Position)
            break
        else
            roi = [roi; temp];
            points = [points; temp.Position];
        end
    end
    oPoints = points;
    
    %draw ROIs
    movAvg = 5;
    slope = 1;
    points = oPoints;
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

%%
    roisignal = zeros(size(img,3),size(points,1));
    for i=1:size(points,1)
        h = drawellipse('Center',oPoints(i,:),'RotationAngle',angles(i),'SemiAxes',[3 3]);
        mask = h.createMask;
        I{i} = find(mask);
    end

    for i=1:size(img,3)
        temp = img(:,:,i);
        for j=1:size(I,2)
            tempI = I{j};
            roisignal(i,j) = mean(temp(tempI),'all');
        end
    end

    roisignalNorm = (roisignal - median(roisignal,1)) ./ median(roisignal,1);
    roismooth = zeros(size(roisignalNorm));
    for i=1:size(roisignalNorm,2)
        roismooth(:,i) = smooth(roisignalNorm(:,i));
    end

    medianroi = median(roismooth,1);
    stdroi = std(roismooth,1);
    thrroi = roismooth > medianroi + 2*stdroi;

    generateActivityMovieIHCs(Kalman_Stack_Filter(single(img)),thrroi,oPoints,angles,[fp '\' name '_activityMovie'],[500 25000])

    IHCstruct = struct();
    IHCstruct.roisignal = roisignal;
    IHCstruct.smoothroi = roismooth;
    IHCstruct.thrroi = thrroi;
    IHCstruct.angles = angles;
    IHCstruct.roiPoints = oPoints;
    save([fp '\' name '_IHCstruct.mat'],'IHCstruct');

end

