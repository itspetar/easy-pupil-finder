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

function[adjXYZNB] = apply_locBias(mleXYZNB,pupilFunc,EMgain,currCoords_rc)
% This function applies localization biases from a pupilFunc structure to
% the MLE fits obtained from localize_emitter_mleInterp.m.

% INPUTS:
% mleXYZNB - 5-element vector of fits before bias correction
% pupilFunc - pupil function structure obtained from ESTIMATE_PUPILFUNC.m
% EMgain - electron multiplication gain used in measurement
% currCoords_rc - coordinates of ROI

% OUTPUTS:
% adjXYZNB - 5-element vector of bias-corrected fits

adjXYZNB = nan(1,5);

adjXYZNB(3) = polyval(pupilFunc.polyDeltaZ,mleXYZNB(3),pupilFunc.errZ,pupilFunc.muZ); % z correction
adjXYZNB(1) = mleXYZNB(1) - polyval(pupilFunc.polyDeltaX,adjXYZNB(3),pupilFunc.errX,pupilFunc.muX)...
    + currCoords_rc(2)*pupilFunc.pixSize; % x correction, based on corrected z
adjXYZNB(2) = mleXYZNB(2) - polyval(pupilFunc.polyDeltaY,adjXYZNB(3),pupilFunc.errY,pupilFunc.muY)...
    + currCoords_rc(1)*pupilFunc.pixSize; % y correction, based on corrected z

adjXYZNB(4) = mleXYZNB(4)*(1 + (EMgain > 1))/polyval(pupilFunc.polyFracN,mleXYZNB(3),pupilFunc.errN,pupilFunc.muN); % N correction
adjXYZNB(5) = mleXYZNB(5)*(1 + (EMgain > 1)); % b correction

end