% Copyright (c)2017, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names 
% of its contributors may be used to endorse or promote products derived 
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function[coords_r,coords_c] = choose_roi(frame,boxSize,title_str,fixedAspectRatioYN)
% This function prompts the user to identify a box within an image. The
% output is the set of coordinates of the rows and columns of the image
% that are chosen.

% INPUTS:
% frame - image data (2D array)
% boxSize - integer corresponding to initial box size in pixels
% title_str - string to be displayed in the figure title
% fixedAspectRatioYN - set to 1 if the box must be a square; otherwise set to 0

% OUTPUTS:
% coords_r - coordinates of rows within 'frame' which were selected
% coords_c - coordinates of columns within 'frame' which were selected

% frame(coords_r,coords_c) yields a 2D array corresponding to the selection
% made within this function.

screenSize = get(0,'Screensize');

selectionFig = figure('Position',screenSize);

subplot(1,2,1)
imagesc(frame)
suptitle(title_str)
colormap(hot);axis image off

selectionBox = imrect(gca,[1,1,boxSize,boxSize]);
setFixedAspectRatioMode(selectionBox,fixedAspectRatioYN) % fix aspect ratio? (1=yes)
setPositionConstraintFcn(selectionBox,makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim')));

badCrop = 1;
while badCrop
    boxLoc = wait(selectionBox);
    boxLoc = round(boxLoc);
    coords_r = boxLoc(2):boxLoc(2)+boxLoc(4)-1;
    coords_c = boxLoc(1):boxLoc(1)+boxLoc(3)-1;
    
    subplot(1,2,2)
    cla
    imagesc(frame(coords_r,coords_c));
    hold on
    plot(round(length(coords_c)/2),round(length(coords_r)/2),'g*')
    colormap(hot);axis image off;title(['zoomed on selection: r=[' num2str(coords_r(1)) ':' num2str(coords_r(end))...
        '], c=[' num2str(coords_c(1)) ':' num2str(coords_c(end)) ']']);
    drawnow
    
    cropPrompt = questdlg('Are you happy with the cropped region?','Cropping success?','Yes','No','Yes');
    switch cropPrompt
        case 'Yes'
            badCrop = 0;
            close(selectionFig)
        case 'No'
            imagesc(zeros(1));axis image off;title('no current selection');
            clear boxLoc
            drawnow
    end
end
end