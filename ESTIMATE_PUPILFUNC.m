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

function[pupilFunc] = ESTIMATE_PUPILFUNC()
% This function carries out a phase retrieval algorithm for engineered point
% spread functions (PSFs) based on maximum likelihood estimation (MLE) of a
% phase aberration term which is added to the theoretical pupil function of
% an imaging system.

% OUTPUT:
% pupilFunc - structure containing pupil function estimation results

screenSize = get(0,'Screensize');

%% 1 - GET FILES AND LOCATIONS FROM USER
[fileName,dirName] = uigetfile('*.tif','Choose stack for MLE phase retrieval');
[darkFileName,darkDirName] = uigetfile('*.tif','Choose dark counts for phase retrieval stack (choose cancel to skip)');
[maskFileName,maskDirName] = uigetfile('*.mat','Load the un-aberrated mask');
saveDir = uigetdir(pwd,'Choose directory in which to save results');

%% 2 - GET MASK PARAMETERS
disp('Preparing phase mask...')
load([maskDirName maskFileName]) % physMask [m]
pupilFunc.physMask = physMask;
maskPrompt1 = questdlg('Is this mask transmissive or reflective?',...
    'Mask file type','Transmissive [m]','Reflective [m]','Transmissive [m]');
switch maskPrompt1
    case 'Transmissive [m]'
        maskPrompt2 = inputdlg({
            'mask array scale [px/m]','refractive index of dielectric mask'...
            },'Provide dielectric mask parameters',1,{'1e5','1.4585'});
        pupilFunc.maskType = 'transmissive';
        pupilFunc.n_mask = str2double(maskPrompt2{2});
    case 'Reflective [m]'
        maskPrompt2 = inputdlg({
            'mask array scale [px/m]'},'Provide phase mask parameters',1,{'1e5'});
        pupilFunc.maskType = 'reflective';
        pupilFunc.n_mask = -1; % effective refractive index for reflective mask (assuming ~90 degree incidence)
end
pupilFunc.px_per_m = str2double(maskPrompt2{1});
disp('... Phase mask ready.')

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

%% 4 - CHOOSE EMITTER AND REPRESENTATIVE BACKGROUND REGION
disp('Identifying signal and background regions...')
frame1 = double(imread([dirName,fileName],1)) - meanDarkImg;
flatStack = frame1;
for ff = 2:nFrames
    flatStack = flatStack + double(imread([dirName,fileName],ff)) - meanDarkImg; % all frames are summed together here
end

if size(frame1,1)>=64
    boxSize = 64;
else
    boxSize = size(frame1,1);
end

[cropCoords_r,cropCoords_c] = choose_roi(frame1,boxSize,'Center box on emitter PSF, then double-click',1);
disp('... Signal region identified.')
[cropCoords_rBG,cropCoords_cBG] = choose_roi(flatStack,boxSize,'Choose representative region for background estimation, then double-click',0);
disp('... Background region identified.')

%% 5 - GET EXPERIMENTAL PARAMETERS FROM USER, IMPORT DATA
disp('Getting experimental parameters...')
experimentPrompt = {
    'initial z0 (emitter position relative to focal plane) [m]',...
    'z step size [m]',...
    'final  z0 (emitter position relative to focal plane) [m]',...
    'measurement EM gain',...
    'emission wavelength [m]',...
    'imaging medium refractive index',...
    'objective lens NA',...
    'microscope magnification',...
    '4f lens focal length [m]',...
    'pixel size [m]',...
    'camera conversion gain'
    };
defaultExperimentInputs = {'3.5e-6','-250e-9','-3.5e-6','100','610e-9','1.518','1.4','133.333','120e-3','120e-9','26.7'};
experimentInputs = inputdlg(experimentPrompt,'Input experimental parameters',1,defaultExperimentInputs);

zStart = str2double(experimentInputs{1});
zStep = str2double(experimentInputs{2});
zEnd = str2double(experimentInputs{3});
EMgain = str2double(experimentInputs{4});

pupilFunc.lambda = str2double(experimentInputs{5});
pupilFunc.n = str2double(experimentInputs{6});
pupilFunc.NA = str2double(experimentInputs{7});
pupilFunc.M = str2double(experimentInputs{8});
pupilFunc.f_4f = str2double(experimentInputs{9});
pupilFunc.pixSize = str2double(experimentInputs{10});
pupilFunc.convGain = str2double(experimentInputs{11});
pupilFunc.zVec = zStart:zStep:zEnd; % vector of z positions

disp('... Experimental parameters imported.')

if length(pupilFunc.zVec) ~= nFrames
    error('Number of axial positions not matched correctly to number of frames of data.')
elseif EMgain < 1
    error('EM gain must be greater than or equal to 1.')
end

% calculate the phase pattern from the physical mask, wavelength, and mask refractive index
pupilFunc.phaseMask = calc_phase_fromPhysMask(pupilFunc.physMask,pupilFunc.lambda,pupilFunc.n_mask);

% import data
disp('Importing data...')
expData = nan([length(cropCoords_r),length(cropCoords_c),nFrames]);
for ii = 1:nFrames
    currFrame = (double(imread([dirName fileName],ii)) - meanDarkImg)*pupilFunc.convGain/EMgain; % convert from counts to photons
    expData(:,:,ii) = currFrame(cropCoords_r,cropCoords_c);
end
expData = max(expData,0); % set negative photon counts to zero
disp('... Data imported.')

%% 6 - ESTIMATE BACKGROUND AND SIGNAL
disp('Calculating background...')
bgData = frame1(cropCoords_rBG,cropCoords_cBG)*pupilFunc.convGain/EMgain; % convert from counts to photons
bgEst = mean(bgData(:)); % calculate mean
bgVec = bgEst*ones(1,nFrames); % [photons]
disp('... Background calculated.')

disp('Calculating signal...')
sigEst = sum(reshape(expData - bgEst*ones(size(expData)),[],nFrames)); % calculate total signal photons in each slice
sVec = sigEst(pupilFunc.zVec == min(abs(pupilFunc.zVec)))*ones(1,nFrames); % [photons], use PSF nearest to z=0
disp('... Signal calculated.')

%% 7 - GET ALGORITHM FEATURES FROM USER
disp('Getting algorithm features...')
algorithmPrompt = {
    'mask angle [deg]',...
    'number of orders of Zernike polynomials',...
    'allow lateral misalignment of mask? [yes/no]'
    };
defaultAlgorithmInputs = {'unknown','4','yes'};
algorithmInputs = inputdlg(algorithmPrompt,'Choose algorithm features',1,defaultAlgorithmInputs);

if strcmp(algorithmInputs{1},'unknown')
    angleImg = expData(:,:,pupilFunc.zVec == max(pupilFunc.zVec));
    pupilFunc.maskAngle = calc_tetrapod_angle(angleImg);
else
    pupilFunc.maskAngle = str2double(algorithmInputs{1});
    if isnan(pupilFunc.maskAngle)
        error('Enter a number for mask angle, or enter ''unknown'' to determine tetrapod mask angle.')
    end
end
pupilFunc.zernikeOrders = str2double(algorithmInputs{2});
nZernikeCoeffs = (pupilFunc.zernikeOrders + 1)*(pupilFunc.zernikeOrders + 2)/2;
pupilFunc.coeffVecStart = [2*rand(nZernikeCoeffs,1) - 1;1;1]; % last two coefficients are lateral mask displacement (x and y, respectively) in pixels
pupilFunc.coeffVecStart(abs(pupilFunc.coeffVecStart) < 0.1) = 0.1; % initial coefficient values are at least 0.1 for all coefficients to avoid problems with fmincon
if strcmp(algorithmInputs{3},'yes') % shiftFlag indicates whether lateral shift of mask is permitted in this optimization
    pupilFunc.shiftFlag = 1;
else
    pupilFunc.shiftFlag = 0;
end

% features which can be tuned if needed:
gBlur = 0.5; % sigma of Gaussian blur filter (see end of get_img_fromCoeffVec)
resizeFactor = 1/4; % reciprocal of oversampling factor in image plane

disp('... Algorithm features imported.')

%% 8 - CALCULATE COORDINATE SYSTEMS
disp('Calculating coordinate systems...')
pupilFunc.FPrad_m = pupilFunc.f_4f*pupilFunc.NA/sqrt(pupilFunc.M^2 - pupilFunc.NA^2); % radius of E-field in Fourier plane (region of E-field support) [m]

[xiGrid_px,etaGrid_px] = meshgrid((1:size(pupilFunc.physMask,1)) - ceil(size(pupilFunc.physMask,1)/2),...
    (1:size(pupilFunc.physMask,1)) - ceil(size(pupilFunc.physMask,1)/2)); % Cartesian coordinates in Fourier plane in units of px, with center (0,0)
[phi,rho] = cart2pol(xiGrid_px/pupilFunc.px_per_m/pupilFunc.FPrad_m,...
    etaGrid_px/pupilFunc.px_per_m/pupilFunc.FPrad_m); % polar coordinates in Fourier plane with rho=1 at edge of E-field support

disp('... Coordinate systems calculated.')

% store coordinate system in pupilFunc structure
pupilFunc.xiGrid_px = xiGrid_px;
pupilFunc.etaGrid_px = etaGrid_px;
pupilFunc.phi = phi;
pupilFunc.rho = rho;

%% 9 - REGISTER X-Y
FOVd = size(expData,1); % image size [px]

% register using unaberrated model
disp('Registering data and imaging model...')
x0y0Vec = nan(nFrames,2); 
regFig = figure;
for zInd = 1:nFrames
    [initModel,~,~,~] = calc_img_fromCoeffVec(zeros(size(pupilFunc.coeffVecStart)),[0,0],pupilFunc.zVec(zInd),...
        1e4,0,FOVd,resizeFactor,gBlur,pupilFunc);
    
    registrationOutput = dftregistration(fft2(initModel),fft2(expData(:,:,zInd)),30);
    rShift = registrationOutput(3);
    cShift = registrationOutput(4);
    x0y0 = [-cShift,-rShift]*pupilFunc.pixSize;
    
    % plot data and shifted (unaberrated) model
    [shiftedModel,~,~,~] = calc_img_fromCoeffVec(zeros(size(pupilFunc.coeffVecStart)),x0y0,pupilFunc.zVec(zInd),...
        1e4,0,FOVd,resizeFactor,gBlur,pupilFunc);
    
    set(0,'CurrentFigure',regFig)
    subplot(1,2,1);
    imagesc(expData(:,:,zInd));
    colormap(hot);axis square off;title('raw data');
    subplot(1,2,2);
    imagesc(shiftedModel);
    colormap(hot);axis square off;title(['model z=' num2str(pupilFunc.zVec(zInd))]);
    pause(0.1)
    
    x0y0Vec(zInd,:) = x0y0;
end
close(regFig)
mean_x0y0 = mean(x0y0Vec); % mean (x,y) shift applied as a global lateral offset
disp('... Registration completed.')

%% 10 - NORMALIZE ALL FIT PARAMETERS
disp('Normalizing fit parameters...')
coeffFactor = 100; % value of 100 represents coefficient value of 1
coeffVecStart_norm = pupilFunc.coeffVecStart*coeffFactor; % normalized version of coeffVecStart
disp('... Parameters normalized.')

%% 11 - SET OPTIMIZATION OPTIONS
options = optimset('MaxFunEvals',1000,'TolFun',1e-1,'TolCon',1,'InitBarrierParam',...
    100,'DiffMinChange',1,'Display','iter','Algorithm','interior-point');

%% 12 - RUN OPTIMIZATION
disp('Beginning optimization...')

% apply excess noise factor if gain is on (in which case optData is in photons/2)
optData = expData./(1 + (EMgain > 1));
optSig = sVec./(1 + (EMgain > 1));
optBG = bgVec./(1 + (EMgain > 1));

% negative log likelihood
llfun = @(x) calc_nLogLikelihood(optData,calc_img_fromCoeffVec(x./coeffFactor,mean_x0y0,pupilFunc.zVec,...
    optSig,optBG,FOVd,resizeFactor,gBlur,pupilFunc));

% constrained minimization
lowerLimVec = [-500*ones(length(pupilFunc.coeffVecStart)-2,1);-3000;-3000];
upperLimVec = [500*ones(length(pupilFunc.coeffVecStart)-2,1);3000;3000];
[resVec,objVal] = fmincon(llfun,coeffVecStart_norm,[],[],[],[],lowerLimVec,upperLimVec,[],options);

% catch errors
if isempty(objVal)
    disp('WARNING: objVal is empty')
    objVal = 0;
end
if or(any(abs(resVec(:) - lowerLimVec(:)) < 0.1),any(abs(resVec(:) - upperLimVec(:)) < 0.1))
    disp('WARNING: optimization hit a boundary.')
end

% get coeffVec_final in coefficient units and radians
coeffVec_final = resVec./coeffFactor;
load('zernike_rad_per_coeff_14orders.mat')
coeffVec_final_rad = coeffVec_final(1:end-2).*transpose(deltaVec(1:nZernikeCoeffs));

% resulting image cube and pupil function (this cube is in photons/2 and has zero background!)
[imgsOut,pupilField,aberrPhase,maskOut] = calc_img_fromCoeffVec(coeffVec_final,mean_x0y0,pupilFunc.zVec,...
    optSig,zeros(size(optBG)),FOVd,resizeFactor,gBlur,pupilFunc);
disp('Optimization completed.')

%% 13 - UPDATE PUPILFUNC STRUCTURE
disp('Storing pupil function results...')
pupilFunc.fileLoc = [dirName,fileName];
pupilFunc.darkFileLoc = [darkDirName,darkFileName];
pupilFunc.maskFileLoc = [maskDirName,maskFileName];
pupilFunc.cropCoords_r = cropCoords_r;
pupilFunc.cropCoords_c = cropCoords_c;
pupilFunc.cropCoords_rBG = cropCoords_rBG;
pupilFunc.cropCoords_cBG = cropCoords_cBG;
pupilFunc.pupilField = pupilField; % full Fourier plane electric field
pupilFunc.aberrPhase = aberrPhase; % aberration term only
pupilFunc.maskOut = maskOut; % mask only (rotated)
pupilFunc.coeffVec = coeffVec_final; % final aberration coefficients (c_j); if included, lateral shifts are last 2 terms
pupilFunc.coeffVec_rad = coeffVec_final_rad; % peak-to-valley aberration size in radians for each mode
pupilFunc.objVal_pfOpt = objVal; % pupil function optimization objective function value
pupilFunc.x0y0Vec = x0y0Vec;
pupilFunc.gBlur = gBlur;
pupilFunc.resizeFactor = resizeFactor;
disp('... Results stored.')

%% 14 - DETERMINE RESIDUAL XYZ BIASES
disp('Determining residual x,y,z biases...')
biasFits = nan(nFrames,5);

% fit PSF data with phase-retrieved model
viewFig = figure('Position',[screenSize(3:4)./4,screenSize(3:4)./2]);
for frameNo = 1:nFrames
    % optData is in units of photons/2, as is optData
    [biasResVec,optImg,~,~,~,~,~] = localize_emitter_mleInterp(optData(:,:,frameNo),imgsOut,pupilFunc.zVec,pupilFunc.zVec(frameNo),pupilFunc.pixSize);
    biasFits(frameNo,:) = biasResVec; % (x,y,z,N,b)

    subplot(1,2,1)
    imagesc(optData(:,:,frameNo))
    axis square off
    colormap hot
    title('raw data')
    subplot(1,2,2)
    imagesc(optImg)
    axis square off
    colormap hot
    title('MLE fit')
    suptitle(['Bias fitting in progress... (' num2str(frameNo) '/' num2str(nFrames) ')'])
    drawnow
end
close(viewFig)
disp('... x,y,z biases determined.')

%% 15 - PLOT AND FIT BIAS RESULTS
disp('Plotting bias results...')
biasFig = figure('Position',screenSize);
badFitOrder = 1;

deltaX = biasFits(:,1);
deltaY = biasFits(:,2);

while badFitOrder
    polyfitDlg = inputdlg({'polynomial fit order:'},'Choose polynomial fit order',1,{'9'});
    fitOrder = str2double(polyfitDlg{1});
    zAxis = linspace(pupilFunc.zVec(1),pupilFunc.zVec(end),100);
    
    % polynomial fits
    [polyDeltaX,errorStructX,muX] = polyfit(pupilFunc.zVec',deltaX,fitOrder);
    [polyDeltaY,errorStructY,muY] = polyfit(pupilFunc.zVec',deltaY,fitOrder);
    [polyDeltaZ,errorStructZ,muZ] = polyfit(biasFits(:,3),pupilFunc.zVec',fitOrder);
    
    % plot lateral and axial biases together with polynomial fits
    set(0,'CurrentFigure',biasFig)
    subplot(2,1,1)
    hold on
    plot(pupilFunc.zVec,deltaX,'b*','linewidth',2)
    plot(pupilFunc.zVec,deltaY,'r*','linewidth',2)
    plot(zAxis,polyval(polyDeltaX,zAxis,errorStructX,muX),'b')
    plot(zAxis,polyval(polyDeltaY,zAxis,errorStructY,muY),'r')
    legend('x','y','polyfit x','polyfit y','Location','northeastoutside')
    xlabel('true z [m]')
    ylabel('bias [m]')
    subplot(2,1,2)
    hold on
    plot(biasFits(:,3),pupilFunc.zVec,'g*','linewidth',2)
    plot(zAxis,polyval(polyDeltaZ,zAxis,errorStructZ,muZ),'g')
    legend('z','polyfit z','Location','southwestoutside')
    xlabel('fitted z [m]')
    ylabel('true z [m]')
    
    suptitle('Press any key to continue...')
    pause
    
    polyfitYN = questdlg({'Polynomial fits OK?'},'Polynomial fits OK?','Yes','No','Yes');
    switch polyfitYN
        case 'Yes'
            badFitOrder = 0;
        case 'No'
            badFitOrder = 1;
            savefig(biasFig,[saveDir '\xyz_bias_results.fig'])
            disp('... Bias results figure saved.')
            clf(biasFig)
    end
end
close(biasFig)
disp('... Plotting completed.')

% calibration of signal photons - used to rescale signal estimate artifacts due to clipping of FOV
[sigModel,~,~,~] = calc_img_fromCoeffVec(coeffVec_final,mean_x0y0,pupilFunc.zVec,...
    zeros(size(optSig)),zeros(size(optBG)),FOVd,resizeFactor,gBlur,pupilFunc);
fracN = squeeze(sum(sum(sigModel,1),2));
fracN = fracN/fracN(pupilFunc.zVec == min(abs(pupilFunc.zVec)));
[polyFracN,errorStructN,muN] = polyfit(pupilFunc.zVec',fracN,fitOrder);

%% 16 - UPDATE PUPILFUNC STRUCTURE
disp('Including biases in pupil function...')
pupilFunc.polyfitOrder =  fitOrder;
pupilFunc.polyDeltaX = polyDeltaX;
pupilFunc.errX = errorStructX;
pupilFunc.muX = muX;
pupilFunc.polyDeltaY = polyDeltaY;
pupilFunc.errY = errorStructY;
pupilFunc.muY = muY;
pupilFunc.polyDeltaZ = polyDeltaZ;
pupilFunc.errZ = errorStructZ;
pupilFunc.muZ = muZ;
pupilFunc.polyFracN = polyFracN;
pupilFunc.errN = errorStructN;
pupilFunc.muN = muN;
pupilFunc.zAxis = zAxis;
pupilFunc.biasFits = biasFits;
disp('... Pupil function updated.')

%% 17 - PLOT PHASE RETRIEVAL RESULTS
disp('Plotting phase retrieval results...')

% plot PSFs (data and phase-retrieved model)
dataFig = figure('Position',screenSize);
set(gcf,'color','white')
simFig = figure('Position',screenSize);
set(gcf,'color','white')
for jj = 1:size(imgsOut,3)
    set(0,'CurrentFigure',dataFig)
    subplot(4,ceil(size(imgsOut,3)/4),jj)
    imagesc(expData(:,:,jj)) % expData is in photons
    title(['raw ' num2str(jj)])
    axis square off
    colormap hot
    
    set(0,'CurrentFigure',simFig)
    subplot(4,ceil(size(imgsOut,3)/4),jj)
    imagesc(imgsOut(:,:,jj).*(1 + (EMgain > 1))) % multiply by 2 to convert to photons if EM gain is used
    title(['sim ' num2str(jj)])
    axis square off
    colormap hot
end
savefig(dataFig,[saveDir '\raw_psf_data.fig'])
savefig(simFig,[saveDir '\mlepr_psf_model.fig'])
disp('... PSF image figures saved.')
close(dataFig)
close(simFig)

% plot pupil function and aberration
pupilFig = plot_field(pupilField,'pupil function',screenSize.*[100,50,0.9,0.35]);
aberrFig = plot_field(aberrPhase,'aberration',screenSize.*[100,screenSize(4)/2 + 50,0.9,0.35]);
savefig(pupilFig,[saveDir '\pupil_field.fig'])
savefig(aberrFig,[saveDir '\phase_aberration.fig'])
disp('... Pupil field and aberration phase figures saved.')
close(pupilFig)
close(aberrFig)

% plot aberration decomposition into coefficients
zernikeFig = figure('Position',screenSize.*[100,100,0.5,0.5]);
set(0,'CurrentFigure',zernikeFig)
plot(4:nZernikeCoeffs,pupilFunc.coeffVec_rad(4:end),'k*')
hline = refline([0,0]);
set(hline,'LineStyle','--')
set(hline,'Color',[0.7,0.7,0.7]);
xlim([3.5,nZernikeCoeffs + 0.5])
set(gca,'XTick',4:nZernikeCoeffs)
ylim([min(coeffVec_final_rad)-0.5,max(coeffVec_final_rad)+0.5])
xlabel('Noll index')
ylabel('peak-to-valley phase [rad]')
set(gcf,'color','white')
if pupilFunc.shiftFlag
    title(['aberration decomposition into Zernike modes (lateral shift [' num2str(pupilFunc.coeffVec(end-1)) ',' num2str(pupilFunc.coeffVec(end)) '] px)'])
else
    title('aberration decomposition into Zernike modes')
end
savefig(zernikeFig,[saveDir '\zernike_decomposition.fig'])
disp('... Zernike decomposition figure saved.')
close(zernikeFig)

disp('... Plotting completed.')

%% 18 - SAVE PUPIL FUNCTION
disp('Saving pupil function...')
save([saveDir '\pupilFunc.mat'],'pupilFunc')
disp('...Pupil function saved.')

disp('PHASE RETRIEVAL COMPLETE.')
end

function[maskAngle] = calc_tetrapod_angle(angleImg)
% This function determines the angle at which a tetrapod mask has been
% placed in an optical system based on an image. It calls the helper
% function lobeFit_ellipticalGauss (below).

% INPUTS:
% angleImg - 2D array with image of tetrapod (should have largest possible value)

% OUPUTS:
% maskAngle - rotation angle needed [degrees]

screenSize = get(0,'Screensize');

badFit = 1;
while badFit
    angleFitFig = figure('Position',[100 100 2*screenSize(3:4)/3]);
    set(0,'CurrentFigure',angleFitFig)
    imagesc(angleImg)
    axis square
    colormap hot
    title('Select lobes, then press enter')

    lobeLocs = ginput; % get lobe locations from user clicks

    % fit each lobe to a single elliptical Gaussian (arbitrary angle)
    [xFit1,yFit1,xSigFit1,ySigFit1,angleFit1,ampFit1,bgFit1] = lobeFit_ellipticalGauss(angleImg,lobeLocs(1,1),lobeLocs(1,2));
    [xFit2,yFit2,xSigFit2,ySigFit2,angleFit2,ampFit2,bgFit2] = lobeFit_ellipticalGauss(angleImg,lobeLocs(2,1),lobeLocs(2,2));

    % check the fits
    [xGrid,yGrid] = meshgrid(1:size(angleImg,2),1:size(angleImg,1));
    lobe1 = ampFit1/(2*pi*xSigFit1*ySigFit1)*exp(-(cos(angleFit1)^2/(2*xSigFit1^2) + sin(angleFit1)^2/(2*ySigFit1^2)).*(xGrid-xFit1).^2 -...
        2*(-sin(2*angleFit1)/(4*xSigFit1^2) + sin(2*angleFit1)/(4*ySigFit1^2)).*(xGrid-xFit1).*(yGrid-yFit1) -...
        (sin(angleFit1)^2/(2*xSigFit1^2) + cos(angleFit1)^2/(2*ySigFit1^2)).*(yGrid-yFit1).^2);
    lobe2 = ampFit2/(2*pi*xSigFit2*ySigFit2)*exp(-(cos(angleFit2)^2/(2*xSigFit2^2) + sin(angleFit2)^2/(2*ySigFit2^2)).*(xGrid-xFit2).^2 -...
        2*(-sin(2*angleFit2)/(4*xSigFit2^2) + sin(2*angleFit2)/(4*ySigFit2^2)).*(xGrid-xFit2).*(yGrid-yFit2) -...
        (sin(angleFit2)^2/(2*xSigFit2^2) + cos(angleFit2)^2/(2*ySigFit2^2)).*(yGrid-yFit2).^2);
    % calculate angle
    maskAngle = rad2deg(atan(abs(yFit2-yFit1)/abs(xFit2-xFit1)));

    hold on
    contour(lobe1 + lobe2,2,'LineColor','blue','LineWidth',3)
    title(['angle = ' num2str(maskAngle)]) % calculated angle appears here
    colormap hot
    axis square

    % decide whether to keep fits
    listFits = {['xFits = ' num2str([xFit1 xFit2])]...
        ['yFits = ' num2str([yFit1 yFit2])]...
        ['xSigFits = ' num2str([xSigFit1 xSigFit2])]...
        ['ySigFits = ' num2str([ySigFit1 ySigFit2])]...
        ['Elliptical angle = ' num2str(rad2deg([angleFit1 angleFit2]))]... % angle of rotation of ellipse w.r.t. x axis
        ['ampFits = ' num2str([ampFit1 ampFit2])]...
        ['bgFits = ' num2str([bgFit1 bgFit2])]...
        ['MASK ANGLE = ' num2str(maskAngle)]};
    fitYN = questdlg(listFits,'Fit OK?','Keep fit','Discard fit','Keep fit');
    switch fitYN
        case 'Keep fit'
            badFit = 0;
        case 'Discard fit'
            badFit = 1;
    end
end

close(angleFitFig)

end

function[xFit,yFit,xSigFit,ySigFit,angleFit,ampFit,bgFit] = lobeFit_ellipticalGauss(img,x0,y0)
% This function is used by calc_tetrapod_angle to fit the lobes of a
% tetrapod PSF to elliptical Gaussians.

init = [x0,y0,4,4,pi,max(img(:)),100];
lb = [x0-3,y0-3,0.25,0.25,0,max(img(:))/3,0];
ub = [x0+3,y0+3,10,10,2*pi,max(img(:))*2,1e4];

[xGrid,yGrid] = meshgrid(1:size(img,2),1:size(img,1));

p = lsqnonlin(@(p) p(6)/(2*pi*p(3)*p(4))*exp(-(cos(p(5))^2/(2*p(3)^2) + sin(p(5))^2/(2*p(4)^2)).*(xGrid-p(1)).^2 -...
    2.*(-sin(2*p(5))/(4*p(3)^2) + sin(2*p(5))/(4*p(4)^2)).*(xGrid-p(1)).*(yGrid-p(2)) -...
    (sin(p(5))^2/(2*p(3)^2) + cos(p(5))^2/(2*p(4)^2)).*(yGrid-p(2)).^2) + p(7) - img,init,lb,ub);

xFit = p(1);
yFit = p(2);
xSigFit = p(3);
ySigFit = p(4);
angleFit = p(5);
ampFit = p(6);
bgFit = p(7);
end

function [phaseMask] = calc_phase_fromPhysMask(physMask,lambda,n)
% This function converts a mask pattern from depth to phase given one
% operating wavelength lambda. The refractive index is used to calculate
% phase delay. Mask is assumed to be placed in air (for which n = 1).

% INPUTS:
% physMask - 2D array containing relative thickness of material with index n (in meters)
% lambda - wavelength of light passing through mask (in meters)
% n - refractive index of mask material

% OUTPUTS:
% phaseMask - 2D array containing relative phases in mask pattern (in radians)

maskPhase = exp(1i*2*pi*((n - 1).*physMask)/lambda);
maskAmp = ones(size(physMask));

phaseMask = maskAmp.*maskPhase;
end