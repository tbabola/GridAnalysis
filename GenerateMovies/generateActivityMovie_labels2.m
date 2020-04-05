function generateActivityMovie_labels2(img,labels,positiveIndices,name,clims)
%generateActivityMovie This function generates a movie of the activity with the grid.
% Variables: img: movie stack (not dFoF, but bleach corrected)
                  %binarywhole = all ROIs represented as binary (active =
                  %1, non-active = 0
%                 postiveIndices = indices for all ROIs inside outlined
%                 region
%                 name = name of file to output
    if nargin < 5
            clims = [500 5000];
    end
    cmap = colormap('prism');
    figure; imagesc(zeros(800,800)); colormap(gray); truesize;
    for i = 1:size(img,3)
        indx_green = find(labels(:,i)>=1);
        imagesc(img(:,:,i)); caxis(clims); hold on;
        %for j=1:size(positiveIndices,1)
        %plot(positiveIndices(j,1:2:end),positiveIndices(j,2:2:end),'Color','k');
        %end
        if size(indx_green,1) > 0
            for j = 1:size(indx_green,1)
                 colors = cmap(mod(labels(indx_green(j),i),256),:);
                 plot(positiveIndices(indx_green(j),1:2:end),positiveIndices(indx_green(j),2:2:end),'Color',colors);
            end
        end
            M(i) = getframe;
        hold off;
    end
    
    if ~isempty(name)  %write video
        v = VideoWriter([name '.mp4'],'MPEG-4');
        v.Quality = 100;

        v.FrameRate = 10;
        open(v);

        for i=1:size(M,2)
            writeVideo(v,M(i).cdata);
        end
        close(v);
    end
end

