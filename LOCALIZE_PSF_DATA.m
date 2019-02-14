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

function[locResults] = LOCALIZE_PSF_DATA()
% This function uses a pupil function determined via ESTIMATE_PUPILFUNC.m
% to perform localizations within a data set. The regions of interest
% containing PSFs are identified manually by the user.

% OUTPUTS:
% locResults - structure containing important results from the analysis

screenSize = get(0,'Screensize');

%% 1 - GET FILES AND LOCATIONS FROM USER
[fileName,dirName] = uigetfile('*.tif','Choose stack to localize');
[darkFileName,darkDirName] = uigetfile('*.tif','Choose dark counts (choose cancel to skip)');
[pfFileName,pfDirName] = uigetfile('*.mat','Choose pupil function');
saveDir = uigetdir(pwd,'Choose directory in which to save results');

%% 2 - GET EXPERIMENTAL PARAMETERS FROM USER
load([pfDirName,pfFileName]) % load pupilFunc structure

disp('Getting experimental parameters...')
experimentPrompt = {
    'maximum z for model [m]',...
    'z step size [m]',...
    'minimum z for model [m]',...
    'measurement EM gain',...
    'desired ROI size [px]'};
defaultExperimentInputs = {num2str(max(pupilFunc.zVec)),'-100e-9',num2str(min(pupilFunc.zVec)),'100','51'};
experimentInputs = inputdlg(experimentPrompt,'Input experimental parameters',1,defaultExperimentInputs);

zStart = str2double(experimentInputs{1});
zStep = str2double(experimentInputs{2});
if zStep > 0
    zStep = -zStep;
end
zEnd = str2double(experimentInputs{3});
zVec = zStart:zStep:zEnd; % vector of z positions

EMgain = str2double(experimentInputs{4});
boxHalf = floor(str2double(experimentInputs{5})/2);

resizeFactor = 1/4;
gBlur = 0.5;

%% 3 - PREPARE DARK COUNTS FOR SUBTRACTION FROM DATA
disp('Preparing dark counts subtraction...')
stackInfo = imfinfo([dirName fileName]);
hh = stackInfo(1).Height;
ww = stackInfo(1).Width;
nFrames = size(stackInfo,1);

% determine mean dark counts image
if darkFileName ~= 0
    darkInfo = imfinfo([darkDirName darkFileName]);
    if or(darkInfo(1).Width ~= ww,darkInfo(1).Height ~= hh)
        error('Dark stack dimensions are different from data dimensions.')
    else
        meanDarkImg = calc_meanDarkCounts(darkDirName,darkFileName);
        disp('... Mean dark counts calculated.')
    end
else
    meanDarkImg = zeros([hh ww]);
    disp(' ... No dark counts file identified. No offset applied to data.')
end

%% 4 - PREPARE PSF DICTIONARY FROM PUPIL FUNCTION
[psfDic,~,~,~] = calc_img_fromCoeffVec(pupilFunc.coeffVec,[0,0],zVec,1e4*ones(size(zVec)),zeros(size(zVec)),2*boxHalf+1,resizeFactor,gBlur,pupilFunc);

%% 5 - PERFORM LOCALIZATION ANALYSIS
disp('Beginning localization analysis...')

% create figure with which user will interact
frameFig = figure('Position',screenSize.*[screenSize(3)*0.1,screenSize(4)*0.1,0.8,0.8]);

% create cell with fit data
locResults.fits = cell(nFrames,1);

% go through each frame
for frameNo = 1:nFrames
    disp(['Current frame: ' num2str(frameNo) '/' num2str(nFrames)])
    
    % load frame
    currFrame = (double(imread([dirName fileName],frameNo)) - meanDarkImg)*pupilFunc.convGain/EMgain; % convert from counts to photons
    
    % display frame
    set(0,'CurrentFigure',frameFig)
    win1 = subplot(5,4,1:16);
    imagesc(currFrame)
    axis image off
    colormap hot
    title('Click on all emitters, then press enter (or press enter to skip)')
    suptitle(['frame ' num2str(frameNo) '/' num2str(nFrames)])
    
    disp('... Please identify the emitters...')
    
    psfLocs = ginput;
    if ~isempty(psfLocs)
        disp(['... ...' num2str(size(psfLocs,1)) ' emitters identified.'])
        
        locResults.fits{frameNo} = nan(size(psfLocs,1),9); % prepare table for PSF data
        
        for psfNo = 1:size(psfLocs,1)
            disp(['... Analyzing emitter ' num2str(psfNo) '/' num2str(size(psfLocs,1))])
            currLoc = psfLocs(psfNo,:);
            currCoords_r = round((currLoc(2) - boxHalf):(currLoc(2) + boxHalf));
            currCoords_c = round((currLoc(1) - boxHalf):(currLoc(1) + boxHalf));
            
            if any([min(currCoords_r) < 1,...
                    max(currCoords_r) > size(currFrame,1),...
                    min(currCoords_c) < 1,...
                    max(currCoords_c) > size(currFrame,2)])
                disp('... ... Emitter too close to edge of frame (skipped).')
            else
                disp('... ... Fitting...')
                
                currROI = currFrame(currCoords_r,currCoords_c); % pick PSF out of frame
                
                % update plot
                set(0,'CurrentFigure',frameFig)
                axes(win1);
                imagesc(currFrame)
                axis image off
                colormap hot
                hold on
                rectangle('Position',[currCoords_c(1),currCoords_r(1),2*boxHalf+1,2*boxHalf+1],'EdgeColor',[0 1 0],'LineWidth',1)
                hold off
                win2 = subplot(5,4,17);
                imagesc(currROI)
                axis image off
                colormap hot
                title('current ROI')
                drawnow
                
                % coarse fit and update plot
                [bestMatch_ind,maxCorrs_all] = find_dicMatch(currROI,psfDic);
                
                set(0,'CurrentFigure',frameFig)
                win3 = subplot(5,4,18);
                plot(zVec,maxCorrs_all,'r')
                hold on
                plot(zVec(bestMatch_ind),maxCorrs_all(bestMatch_ind),'*')
                xlabel('z match [m]')
                xlim([min(zVec),max(zVec)])
                ylabel('correlation')
                title('correlation with dictionary')
                win4 = subplot(5,4,19);
                imagesc(psfDic(:,:,bestMatch_ind))
                axis image off
                colormap hot
                title('best dictionary match')
                drawnow
                
                % fine fit and apply bias correction
                [mleXYZNB,optImg,~,~,boundFlag,~,bgQuant] = localize_emitter_mleInterp(currROI./(1 + (EMgain > 1)),...
                    psfDic,zVec,zVec(bestMatch_ind),pupilFunc.pixSize);
                adjXYZNB = apply_locBias(mleXYZNB,pupilFunc,EMgain,[currCoords_r(1),currCoords_c(1)]);
                
                % update plot
                set(0,'CurrentFigure',frameFig)
                win5 = subplot(5,4,20);
                imagesc(optImg.*(1 + (EMgain > 1)))
                axis image off
                colormap hot
                title('fine MLE fit')
                drawnow
                
                % store data
                locResults.fits{frameNo}(psfNo,:) = [frameNo,psfNo,adjXYZNB,boundFlag,bgQuant];
                
                % clear plot
                set(0,'CurrentFigure',frameFig)
                cla(win2);cla(win3);cla(win4);cla(win5);
            end
        end
    else
        disp('... ... 0 emitters identified (skipping frame).')
    end
    disp('... Frame analyzed.')
end

disp('... Analysis complete.')
close(frameFig)

%% 6 - UPDATE OUTPUT STRUCTURE
locResults.fileName = fileName;
locResults.dirName = dirName;
locResults.darkFileName = darkFileName;
locResults.darkDirName = darkDirName;
locResults.pfFileName = pfFileName;
locResults.pfDirName = pfDirName;
locResults.pupilFunc = pupilFunc;
locResults.EMgain = EMgain;
locResults.zVec = zVec;
locResults.resizeFactor = resizeFactor;
locResults.gBlur = gBlur;
locResults.psfDic = psfDic;

%% 7 - SAVE RESULTS
disp('Saving localization results...')
save([saveDir '\locResults.mat'],'locResults')
disp('...Localization results saved.')

disp('LOCALIZATION ANALYSIS COMPLETE.')
end