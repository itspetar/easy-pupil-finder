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

function[imgsOut,pupilField,aberrPhase,maskRot] = calc_img_fromCoeffVec(coeffVec,x0y0,zVec,sVec,bgVec,FOVd,resizeFactor,gBlur,pupilFunc)
% This function generates a set of point spread functions at various depths
% within the imaging medium by adding a phase aberration described by a
% combination of Zernike modes to the theoretical electric field in the
% Fourier plane.

% INPUTS:
% coeffVec - column vector of optimization parameters
%     first terms are Zernike coefficients, last two terms are lateral mask shifts
% x0y0 - 2-element array, lateral offsets from the center of the image [m]
% zVec - vector of z positions of emitters [m]
% sVec - vector of signal [total counts]
% bgVec - vector of background intensity [counts/px]
% FOVd - diameter of output images [px]
% resizeFactor - determines oversampling of image plane (lower --> greater oversampling)
% gBlur - standard deviation of Gaussian blur applied to images [px]
% pupilFunc - pupil function structure

% OUTPUTS:
% imgsOut - calculated image plane intensity
% pupilField - Fourier (pupil) plane electric field
% aberrPhase - aberration phase applied to Fourier plane electric field
% maskRot - phase mask following rotation and/or lateral shift

%% APPLY ROTATION AND SHIFT TO MASK
maskRot = imrotate(pupilFunc.phaseMask,pupilFunc.maskAngle,'crop'); % mask rotation

if pupilFunc.shiftFlag % lateral mask shift (if shiftFlag=1)
    xiShift = coeffVec(end-1);
    etaShift = coeffVec(end);
    angleMask = interp2(pupilFunc.xiGrid_px,pupilFunc.etaGrid_px,angle(maskRot),...
        pupilFunc.xiGrid_px + xiShift,pupilFunc.etaGrid_px + etaShift,'cubic');
    absMask = interp2(pupilFunc.xiGrid_px,pupilFunc.etaGrid_px,abs(maskRot),...
        pupilFunc.xiGrid_px + xiShift,pupilFunc.etaGrid_px + etaShift,'cubic');
    maskRot = absMask.*exp(1i*angleMask);
    maskRot(isnan(maskRot)) = 0;
end

coeffVec = coeffVec(1:end-2); % delete shift terms from coeffVec for rest of this function
coeffVec(1) = 0; % zero out piston mode (tip & tilt included for fine lateral alignment during optimization)

%% CALCULATE PUPIL FUNCTION TERMS
ampMod = (1 - (pupilFunc.NA*pupilFunc.rho/pupilFunc.n).^2).^(-1/4); % wavefront compression (apodization)
ampMod(pupilFunc.rho > 1) = 0; % circ function (limits bfp support)
lateralPhase = exp(1i*2*pi*pupilFunc.M*pupilFunc.NA*pupilFunc.rho.*(cos(pupilFunc.phi)*x0y0(1) +...
    sin(pupilFunc.phi)*x0y0(2))/pupilFunc.lambda/sqrt(pupilFunc.M^2 - pupilFunc.NA^2)); % lateral shift

% Zernike phase:
D = calc_zernikeBasis_Noll(round(2*pupilFunc.FPrad_m*pupilFunc.px_per_m),...
    pupilFunc.zernikeOrders(pupilFunc.zernikeOrders > 0),size(maskRot,1));
zernikeAberrs = reshape(D*coeffVec,size(maskRot));
aberrPhase = exp(1i.*zernikeAberrs); % additional phase term

pupilField = ampMod.*aberrPhase.*maskRot.*lateralPhase; % pupil function

%% GENERATE INDIVIDUAL IMAGES ACCORDING TO Z0 VALUES
padVal = round((pupilFunc.lambda*pupilFunc.f_4f*pupilFunc.px_per_m/(resizeFactor*pupilFunc.pixSize*pupilFunc.M)...
    - size(maskRot,1))/2); % fft padding on one side [px]
imgsOut = nan([FOVd FOVd length(zVec)]);

for ii = 1:size(imgsOut,3)
    axialPhase = 2*pi*pupilFunc.n*zVec(ii)*sqrt(1 - (pupilFunc.NA*pupilFunc.rho/pupilFunc.n).^2)/pupilFunc.lambda; % phase due to depth term
    axialPhase(pupilFunc.rho > 1) = 0;
    axialPhase = exp(1i.*axialPhase);
    
    fpField = axialPhase.*pupilField; % complete Fourier plane field term

    % generate image from pupil function:
    Eimg_OS = fftshift(fft2(ifftshift(padarray(fpField,[padVal padVal],'both')))); % oversampled image plane E-field
    Iimg_full = imresize(Eimg_OS.*conj(Eimg_OS),resizeFactor,'bilinear'); % image plane intensity (downsampled to system pixel size)
    Iimg = Iimg_full(round(end/2 - (FOVd - 1)/2):round(end/2 + (FOVd - 1)/2),...
        round(end/2 - (FOVd - 1)/2):round(end/2 + (FOVd - 1)/2)); % crop image
    
    if sum(sVec) ~= 0 % generally true; if signal=0 the images will not be scaled or have background added
        % scale signal, add background:
        Iimg = Iimg*sVec(ii)/sum(Iimg(:)); % rescale signal
        Iimg = Iimg + max(0,bgVec(ii)); % add background
    end

    % add Gaussian blur (optional):
    if gBlur > 0
        h = fspecial('gaussian',[5 5],gBlur);
        Iimg = imfilter(Iimg,h);
    end
    
    % store:
    imgsOut(:,:,ii) = Iimg;
end
end

function [D] = calc_zernikeBasis_Noll(diam,maxOrder,boxLength)
% This function calculates a normalized Zernike basis, arranged according
% to the Noll ordering, as in R. J. Noll, "Zernike polynomials and
% atmospheric turbulence," J. Opt. Soc. Am. 66, 207-211 (1976). The helper
% function get_mInds_Noll (below) is used to calculate the correct sequence
% of polynomials.

% INPUTS:
% diam - diameter of circular support [px]
% maxOrder - highest polynomial order to be included
% boxLength - side length of box [px]

% OUTPUTS:
% D - 2D matrix containing Zernike basis data.

xx = linspace(-1,1,diam);
yy = xx;
[XX,YY] = meshgrid(xx,yy);
rho = sqrt(XX.^2+YY.^2);
phi = atan2(YY,XX);

ind = 0;

for n = 0:maxOrder
    mVec = get_mInds_Noll(n);
    for m = mVec
        ind = ind+1;

        % get radial part
        Rmn = 0;
        absm = abs(m);
        for k=0:(n-absm)/2
            Rmn = Rmn + ((-1)^k)*(factorial(n-k))*rho.^(n-2*k)/((factorial(k))*(factorial((n+absm)/2-k))*(factorial((n-absm)/2-k)));
        end

        % add azimuthal part
        if m > 0
            Zmn = Rmn.*cos(absm.*phi);
        elseif m < 0
            Zmn = Rmn.*sin(absm.*phi);
        else
            Zmn = Rmn;
        end

        Zmn((XX.^2+YY.^2)>1)=0;
        if mod(diam,2) ~= 0
            Zmn = padarray(Zmn,ceil([(boxLength-diam)/2,(boxLength-diam)/2]),'pre');
            Zmn = padarray(Zmn,floor([(boxLength-diam)/2,(boxLength-diam)/2]),'post');
        else
            Zmn = padarray(Zmn,[boxLength-diam,boxLength-diam]/2,'both');
        end

        D(:,ind) = Zmn(:);
    end
end

D = normc(D)*sqrt((diam/2)^2*pi); % normalize to pi

end

function[vec] = get_mInds_Noll(n)
% This function determines all m indices for the nth order of Zernike
% polynomials arranged with Noll ordering (helper function to
% calc_zernikeBasis_Noll).
    prevInd = sum(0:n);

    vec = -fliplr(-n:2:0);
    vec = [vec; -vec]; % [pos; neg]

    if mod(mod(prevInd,2) + mod(n,2),2) == 1
        vec = flipud(vec);
    end

    if mod(n,2) == 0
        vec = transpose(vec(2:end));
    else
        vec = vec(:);
    end

    vec = transpose(vec);
end