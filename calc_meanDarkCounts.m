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

function[meanDarkImg] = calc_meanDarkCounts(darkDirName,darkFileName)
% This function determines the mean dark counts image for subtraction from
% a data set.

% INPUTS:
% darkDirName - path to TIF containing dark counts (string)
% darkFileName - file name of TIF containing dark counts (string)

% OUTPUTS:
% meanDarkImg - mean dark image in units of counts

darkInfo = imfinfo([darkDirName darkFileName]);
nDarkFrames = size(darkInfo,1);
meanDarkImg = zeros(darkInfo(1).Height,darkInfo(1).Width);

for frameNo = 1:nDarkFrames
    dFrame = double(imread([darkDirName darkFileName],frameNo));
    meanDarkImg = meanDarkImg + dFrame;
end
meanDarkImg = meanDarkImg./nDarkFrames;

end