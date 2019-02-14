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

function [h] = plot_field(data,titleStr,posVec)
% This function plots a complex field (such as a pupil function) by
% separating it into an amplitude component and a phase component. It loads
% 'phasemap.mat', which is the phase colormap.

% INPUTS:
% data - complex field data
% titleStr - desired title string in image
% posVec - position vector for figure placement

% OUPUTS:
% h - figure handle

load('phasemap.mat')

if exist('h','var')
    h(end+1) = figure('position',posVec);
else
    h(1) = figure('position',posVec);
end

subplot(1,7,4:7)
imagesc(angle(data))
title('phase')
axis square
colormap(phasemap);
caxis([-pi pi])
cbar1 = colorbar;
set(cbar1,'YTick',[-pi 0 pi])
set(cbar1,'TickLabels',{'-\pi' '0' '\pi'})

subplot(1,7,1:3)
imagesc(abs(data))
title('amplitude')
axis square
colormap(imgca,gray)
colorbar

suptitle(titleStr)
set(gcf,'color','white')
end