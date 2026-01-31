% Sua Bae
%
%  2019-05-22
%  IQ/RF beamformer for plane wave imaging
%  linear
%
function stBFpm = SetBFParam(sBfType, stRFInfo, nMultiBeam, nQiq, nQph, sBfOutName)
    
    stTrans = stRFInfo.stTrans;
    assert(strcmp(stRFInfo.sVersion,'2.0')||strcmp(stRFInfo.sVersion,'2.1'),'version mismatch'); % nRcvChannel should be included
    
    
    stBFpm.sVersion = '3.0';
    stBFpm.bPlot = 0;
    stBFpm.bSaveEachAngle = 0;    
    
    if exist('sBfOutName','var')
        stBFpm.sBfOutName = sBfOutName;
    else
        stBFpm.sBfOutName = 'stBfData.mat'; % default
    end
    
    %% 1. Demodulation/ Decimation/ Interpolation
    
    if strcmp(sBfType,'RFBF')
        stBFpm.nQiq = nQiq; % 2 for 2f0 delay resolution
        if nQiq == 2
            stBFpm.nDeciFactor = 2;
            stBFpm.nIntpFactor = 1;
        elseif nQiq >= 4
            stBFpm.nDeciFactor = 1;
            stBFpm.nIntpFactor = nQiq/4;
        end
        
        stDemod.bDynamic = false;
          stDemod.nDemodFreq       = stRFInfo.nFc; % demodulation frequency
          stDemod.nDemodLPFFreq = stRFInfo.nFs/4; % just remove replica (do not bandlimit)
          stDemod.nDemodLPFTap = 65;
        stBFpm.stDemod = stDemod;
        
    elseif strcmp(sBfType,'IQBF')
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
          stDemod.nDemodLPFFreq = stRFInfo.nFs/4; % just remove replica (do not bandlimit)
          stDemod.nDemodLPFTap = 65;
        stBFpm.stDemod = stDemod;

        stBFpm.nDeciFactor = 2; 
        stBFpm.nIntpFactor = nQiq/4*stBFpm.nDeciFactor; % originally 4f0 sampling
    end

    assert(stBFpm.nIntpFactor>=1,'stBFpm.nIntpFactor is smaller than 1!!!');
    assert(stBFpm.nIntpFactor == double(uint8(stBFpm.nIntpFactor)),'stBFpm.nIntpFactor must be an integer!!!');
    
    
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
    stBFpm.nSysOffset_sec   = 2e-6;      
         
    
    
end