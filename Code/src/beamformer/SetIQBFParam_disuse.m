% Sua Bae
%
% 2019-05-20
%    update: nBW [%]: ex) nBW=20 for 20% BW
%            stDemod.nDemodLPFFreq = stTrans.nFs/4; % just remove replica (do not bandlimit)
%
% 2019-04-27
%    update stDemod.nDemodLPFFreq, stDemod.nDemodLPFTap
%           use uint8 to make sure integer for nIntpFactor, nDeciFactor
% 2019-04-05
%  stBFpm.sVersion 2.1->3.0
%     allow stRFInfo.sVersion 2.1 
%     limit the bandwidth of original signal : nBW (0.5 for 50% BW)
%     control the tiem/phase delay resolution : nQiq, nQph
%
% 2019-03-14
%  stBFpm.sVersion 2.0->2.1
%    include stBFpm.nMultiBeam
%            stBFpm.sBfOutName
%            stBFpm.bDSC
%
% 2019-01-27
%  stRFInfo version 2.0
%      nRcvChannel is used
%      nTxOffsetDelay_sec -> nSysOffset_sec
%
% 2017-12-21
%  import stTrans from stRFInfo
% 2017-11-28
%   view angle: 30 deg for xz plane, 23.3584 deg for yz plane
% 2017-11-09
%  using clBFGrid3D
%
function stBFpm = SetIQBFParam(stRFInfo, nMultiBeam, nQiq, nQph, sBfOutName)
    
    stTrans = stRFInfo.stTrans;
    assert(strcmp(stRFInfo.sVersion,'2.0')||strcmp(stRFInfo.sVersion,'2.1'),'version mismatch'); % nRcvChannel should be included
    
    
    stBFpm.sVersion = '3.0';
    stBFpm.bPlot = 1;
    stBFpm.bSaveEachAngle = 0;    
    
    if exist('sBfOutName','var')
        stBFpm.sBfOutName = sBfOutName;
    else
        stBFpm.sBfOutName = 'stBfData.mat'; % default
    end
    
    %% 1. Demodulation/ Decimation/ Interpolation (IQ Data!)  
    stBFpm.nQiq = nQiq; % 2 for 2f0 delay resolution
    if nQph ~= 0
        stBFpm.nQph = nQph; % 8 for 2*pi/8 phase delay resolution
    end
%     stDemod.bDynamic = true;
%         stDemod.nDemodFreq_slop         = -0.008*1e3*1e6; % [Hz/m]
%         stDemod.nDemodFreq_intercept    = 2.9784*1e6; % [Hz]
%         stDemod.nDemodLPFFreq    = 1e6;%min(stRFInfo.nFc, stRFInfo.nFs/2/stBFpm.nDeciFactor); % demodulation LPF frequency % "nFc" (200% band limit) and "nFs/2/nDeciFactor" (band limit for decimation)
%         stDemod.nDemodLPFTap     = 127; % demodulation LPF tap
    stDemod.bDynamic = false;
      stDemod.nDemodFreq       = stRFInfo.nFc; % demodulation frequency
%       switch nBW
%           case 0.5;  nFactor2 = 0.63;
%           case 0.75;  nFactor2 = 0.725;
%           case 1;  nFactor2 = 0.83;
%           case 1.25;  nFactor2 = 1.25;
%       end        
%       stDemod.nDemodLPFFreq    = nBW*stDemod.nDemodFreq/2*nFactor2; % b*fc/2*nFactor2 (one-side, nFactor2: for -45dB bandwidth instead of -6 dB)
%       stDemod.nDemodLPFTap     = 30/nBW+1; % demodulation LPF tap
      stDemod.nDemodLPFFreq = stRFInfo.nFs/4; % just remove replica (do not bandlimit)
      stDemod.nDemodLPFTap = 65;
    stBFpm.stDemod = stDemod;
    
    stBFpm.nDeciFactor = 2; 
    stBFpm.nIntpFactor = double(uint8(nQiq/4*stBFpm.nDeciFactor)); % originally 4f0 sampling
%     
%     nOrgFactor = stRFInfo.nFs/stRFInfo.nFc;
%     % reduce sampling rate from nFs to around 2*nBW*nFc and must be a even integer
%     stBFpm.nDeciFactor  = double(uint8(round(stRFInfo.nFs/(2*stRFInfo.nBW/100*stRFInfo.nFc)/2)*2)); % -- 2019.05.20
% %     stBFpm.nDeciFactor  = double(uint8(round(stRFInfo.nFs/(2*nBW*stRFInfo.nFc)/2)*2)); -- 2019.04.27
% %     stBFpm.nDeciFactor  = floor(stRFInfo.nFs/(2*nBW*stRFInfo.nFc)); 
% %     stBFpm.nDeciFactor  = floor(stRFInfo.nFs/(nBW*stRFInfo.nFc)); 
% %     stBFpm.nDeciFactor = 1;
%     assert(stBFpm.nDeciFactor>=1,'stBFpm.nDeciFactor is smaller than 1!!!');
%     
%     % increase sampling rate from decimated frequency to nQiq (integer!!!)
%     stBFpm.nIntpFactor = double(uint8(nQiq/(nOrgFactor/stBFpm.nDeciFactor)));
% %     stBFpm.nIntpFactor = 1;
    assert(stBFpm.nIntpFactor>=1,'stBFpm.nIntpFactor is smaller than 1!!!');
    assert(nQiq == 4/stBFpm.nDeciFactor*stBFpm.nIntpFactor,'stBFpm.nIntpFactor is smaller than 1!!!');
    
    
    %% 2. Beamforming grid
    if exist('nMultiBeam','var')
        stBFpm.nMultiBeam = nMultiBeam;
    else
        stBFpm.nMultiBeam = 1; % default
    end
    %%%%--  Dim, interval, offset
    % num of BF points
    stG.nXdim = stTrans.nNumEle*stBFpm.nMultiBeam;
    stG.nZdim = round(stRFInfo.nSample/stBFpm.nDeciFactor*stBFpm.nIntpFactor);
    % stG.nZdim = round(nReconDepth*2/stRFInfo.nSoundSpeed*stRFInfo.nFs/stBFpm.nDeciFactor*stBFpm.nIntpFactor);
    % stG.nZdim = stRFInfo.nSample;
    % interval of BF points
    stG.dx = stTrans.nPitch_x/stBFpm.nMultiBeam; % [deg]
    stG.dz = 1/(stRFInfo.nFs)*(stRFInfo.nSoundSpeed)/2*stBFpm.nDeciFactor/stBFpm.nIntpFactor; % [m] dr = sample interval of RF data 
    % axes
    stG.aX = (-(stG.nXdim-1)/2:1:(stG.nXdim-1)/2)*stG.dx;
    stG.aZ = (0:1:(stG.nZdim-1))*stG.dz;

    %%%%-- BF grid
    display(['ROI size: X: ' num2str(stG.aX(1)*1e3) '~' num2str(stG.aX(end)*1e3) 'mm, ' ...
                       'Z: ' num2str(stG.aZ(1)*1e3) '~' num2str(stG.aZ(end)*1e3) 'mm']);
    display(['          dx: ' num2str(stG.dx*1e3) 'mm (' num2str(stG.dx/(stRFInfo.nSoundSpeed/stRFInfo.nFc)) ' lambda), ' ...
                       'dz: ' num2str(stG.dz*1e3) 'mm (' num2str(stG.dz/(stRFInfo.nSoundSpeed/stRFInfo.nFc)) ' lambda)']);  
                    
    stBFpm.stG = stG;

   %% 3. DSC grid
    stBFpm.bDSC = false; % when using 'Mid Processor', don't need to DSC here
        
    %% 3. RX
    stBFpm.nRxFnum          = 1.0;
    stBFpm.sRxApodWindow    = 'boxcar';%'tukey50';
    
    %% 4. system offset delay
    stBFpm.nSysOffset_sec   = 0e-6;      
         
    
    
end