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

function[bestMatch_ind,maxCorrs_all] = find_dicMatch(img,dic)
% This function finds the best match from a dictionary to an image.

% INPUTS:
% img - 2D array (image)
% dic - 3D array (dictionary; 3rd dimension is searched)

% OUTPUTS:
% bestMatch_ind - index of best match (i.e. dic(:,:,bestMatch_ind) is the best
% maxCorrs_all - vector with length size(dic,3), gives correlation as a function of the 3rd dictionary dimension

corrMaps_all = nan([size(img,1),size(img,2),size(dic,3)]);
maxCorrs_all = nan(1,size(dic,3));

% fill correlation matrix
for ff = 1:size(dic,3)
    modelImg_try = dic(:,:,ff);
    corr_slice = conv2(modelImg_try/norm(modelImg_try(:)),img/norm(img(:)),'same');
    corrMaps_all(:,:,ff) = corr_slice;
    maxCorrs_all(ff) = max(corr_slice(:));
end

bestMatch_ind = find(maxCorrs_all == max(maxCorrs_all));
end