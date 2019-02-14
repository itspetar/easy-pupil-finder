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

function[fitXYZNB,optImg,initXYZNB,objVal,boundFlag,optImg_noBG,bgQuant] = localize_emitter_mleInterp(expData,dic,dicVals,initZ,pixSize)
% This function fits a PSF using maximum likelihood estimation. The lateral
% position, signal photons, background photons, and axial position of the
% emitter are the parameters of the fit.

% INPUTS:
% expData - PSF image to be localized [2D array]
% dic - dictionary of model PSFs used for fit [3D array; 3rd dimension is z]
% dicVals - vector of z positions of dictionary [meters]
% initZ - initial z guess [meters]
% pixSize - pixel size [meters]

% OUTPUTS:
% fitXYZNB - 5-element vector containing fits (x,y,z,signal,background/px)
% optImg - optimal image (calculated with parameters fitXYZNB)
% initXYZNB - 5-element vector containing initial points
% objVal - value of objective function (negative log likelihood) at optimum
% boundFlag - equals 1 (-1) if an upper (lower) bound is hit, 0 otherwise
% optImg_noBG - optimal image as in optImg but with 0 background photons/px
% bgQuant - quantile of image data used as background (fit may be bad if this is >0.3)

%% OPTIMIZATION OPTIONS
options = optimset('MaxFunEvals',1000,'TolFun',1e-3,'TolCon',1e-3,'InitBarrierParam',...
    100,'DiffMinChange',1e-3,'Display','off','Algorithm','interior-point');

%% INITIALIZATION, XY REGISTRATION
if abs(initZ) < pixSize/10
    initZ = initZ + (pixSize/10)*sign(initZ + eps);
end

refVec = {'x','y','z','N','b'};

FOVd = size(expData,1); % size of image

% register using initial guess
initModel = interp_psfDic(dic,dicVals,[0 0],initZ,2e4,0,pixSize);
regOutput = dftregistration(fft2(initModel),fft2(expData),30);
rShift = regOutput(3);
cShift = regOutput(4);

% initial parameters
initXY = [-cShift,-rShift]*pixSize; % lateral position

% catch zeros to avoid conflict with optimizaton routines
if abs(initXY(1)) < pixSize/10
    initXY(1) = initXY(1) + (pixSize/10)*sign(initXY(1) + eps);
end
if abs(initXY(2)) < pixSize/10
    initXY(2) = initXY(2) + (pixSize/10)*sign(initXY(2) + eps);
end

itersOut = 1;
bgQuant = 0;

while itersOut == 1
    bgQuant = bgQuant + 0.1;
    if bgQuant > 0.1
        disp(['increased background quantile to ' num2str(bgQuant)])
    end
    initB = max(quantile(expData(:),bgQuant),1); % background photons estimated from quantile of data
    initN = max(500,sum(expData(:)) - initB*length(expData(:))); % signal photons
    initXYZNB = [initXY,initZ,initN,initB]; % vector of initial parameters

    %% NORMALIZE ALL FIT PARAMETERS TO 100
    scaleFactor = 100./initXYZNB;
    startXYZNB_norm = initXYZNB.*scaleFactor; % normalized 

    %% RUN OPTIMIZATION
    % negative log likelihood
    llfun = @(x) calc_nLogLikelihood(expData,interp_psfDic(dic,dicVals,[x(1)/scaleFactor(1),x(2)/scaleFactor(2)],...
        x(3)/scaleFactor(3),x(4)/scaleFactor(4),x(5)/scaleFactor(5),pixSize));

    % constrained minimization
    lowerLimVec = [-5000;-5000;-5000;0;0];
    upperLimVec = [500;500;5000;1000;1000];
    [resVec,objVal,~,output] = fmincon(llfun,startXYZNB_norm,[],[],[],[],lowerLimVec,upperLimVec,[],options);
    itersOut = output.iterations;
    if bgQuant > 0.8
        itersOut = -1;
    end
end
    
% catch errors
if isempty(objVal)
    disp('objVal empty')
    objVal = 0;
end
if any(abs(resVec(:) - lowerLimVec(:)) < 0.1)
    disp(['Optimization hit lower bound (' refVec{abs(resVec(:) - lowerLimVec(:)) < 0.1} ')'])
    boundFlag = -1;
elseif any(abs(resVec(:) - upperLimVec(:)) < 0.1)
    disp(['Optimization hit lower bound (' refVec{abs(resVec(:) - upperLimVec(:)) < 0.1} ')'])
    boundFlag = 1;
else
    boundFlag = 0;
end

% get fitXYvNb in real units from resVec
fitXYZNB = resVec./scaleFactor;

%% GET RESULTING IMAGE
optImg = interp_psfDic(dic,dicVals,[fitXYZNB(1),fitXYZNB(2)],fitXYZNB(3),fitXYZNB(4),fitXYZNB(5),pixSize);
optImg_noBG = interp_psfDic(dic,dicVals,[fitXYZNB(1),fitXYZNB(2)],fitXYZNB(3),fitXYZNB(4),0,pixSize);
end

function[img] = interp_psfDic(dic,dicVals,x0y0,z0,sig,bg,pixSize)
% This function uses 3D interpolation to calculate a PSF at (x0,y0,z0) from
% a dictionary with x=0, y=0, and z given by dicVals.

% INPUTS:
% dic - 3D array containing the PSF dictionary (third dimension is axial dimension)
% dicVals - axial positions of dictionary slices
% x0y0 - 2-element vector corresponding to desired lateral position [meters]
% z0 - desired axial position [meters]
% sig - total signal within PSF
% bg - constant background [counts/px]
% pixSize - pixel size [meters]

% OUPUTS:
% img - calculated image plane intensity after shifting to (x0,y0,z0)

[~,zVicinityMid] = min(abs(dicVals - z0));
zVicinity = max(1,zVicinityMid - 2):min(zVicinityMid + 2,size(dicVals,2));
[X,Y,Z] = meshgrid(1:size(dic,2),1:size(dic,1),dicVals(zVicinity));
[Xq,Yq,Zq] = meshgrid((1:size(dic,2)) - x0y0(1)/pixSize,(1:size(dic,1)) - x0y0(2)/pixSize,z0);
img = interp3(X,Y,Z,dic(:,:,zVicinity),Xq,Yq,Zq,'cubic',bg);

% scale signal and add background
img = img*sig/sum(img(:));
img = img + max(0,bg);
end