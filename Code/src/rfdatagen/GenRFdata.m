% Sua Bae
% 2019-07-02
%   Version 2.6
%       add Tx apodization if mApod exists in stTx structure
% 2019-04-26
%   Version 2.5
%       remove input: nImpulseRespPulseCycle
%       update: Gaussian Envelop for bandwidth control!!!
%       add a field in stTrans structure: 'stTrans.nBW_45dB' or 'stTrans.nBW_6dB'
% 2019-04-01
%   Version 2.1
%       using TxAptrGen.m(2019-04-01) instead of TxApodGen.m
% 2019-03-29
%   Version 2.0
%       new input: nImpulseRespPulseCycle (to control bandwidth of round-trip signal)
%       removed input: nSoundSpeed (now it is fixed to 1540 m/s)
%       
% 2017-08-25
%   point spread function
%
% 2025-01-31
%   Modified: replaced padarray with basic MATLAB functions (no Image Processing Toolbox needed)
%
function [stRcvData] = GenRFdata(stTrans, stTx, stScat, bPlot)
 
     sVersion = '2.6';
        
     bShowFieldIISetUp = bPlot;

    %% Field II setup
    addpath('C:\Field_II');
    
    %%
    sStartAcqTime = 'first firing time'; 

   
    %%% Ultrasound parameter
    nAttenFactor 	= 0.5/(1e-2*1e6);   % Attenuation factor[dB/(m*Hz)] % 0.5 dB/(cm*MHz)
    nFc             = stTrans.nFc; % [Hz]
    nFs             = stTrans.nFc * 4; % [Hz]
    nUpSampFactor   = 8;
    nSimulFrequency	= nFs*nUpSampFactor; % [Hz]
    nSoundSpeed     = 1540;
    nWaveLength 	= nSoundSpeed/nFc;    
    dt              = 1/nSimulFrequency; % [sec]
    nTc             = 1/nFc;        % period [sec]
    nTxPulseCycle   = stTx.nTxPulseCycle; % default TxPulseCycle is set to 1
    
    %%% Field
    field_init(0);
    set_field('c',nSoundSpeed);             % Set speed of sound
    set_field('Freq_att',nAttenFactor);     % frequency dependency attenuation in [dB/(m*Hz)] around the center frequency
    set_field('att', nFc*nAttenFactor);     % Freqency independent attenuation[dB/m] 
    set_field('att_f0',nFc);                % Attenuation center frequency[Hz]
    set_field('use_att',1);                 % attenuation on/off
    set_field('fs',nSimulFrequency);        % set the sampling frequency
    set_field('use_rectangles',1);          % use ractangles for apertures

    %%%
    aFocalPointPos = [0, 0, 20]*1e-3; % Focal point (not used but need to be defined)


    %%% Transducer
    nFarFieldDepth_x = 0.01e-3; % depth over nFarFieldDepth_x is assumed to be far field on x-z plane 
    nFarFieldDepth_y = 0.5e-3; % depth over nFarFieldDepth_x is assumed to be far field on y-z plane 

    nMathEleSize_x = sqrt(nFarFieldDepth_x * 4 * nWaveLength); % H < sqrt(4*lambda*z)
    nMathEleSize_y = sqrt(nFarFieldDepth_y * 4 * nWaveLength);

    nSubDivNum_x = ceil(stTrans.nWidth_x / nMathEleSize_x);
    nSubDivNum_y = ceil(stTrans.nWidth_y/ nMathEleSize_y);

    if strcmp(stTrans.sType,'Linear')
        % Has a elevationally-focused lens 
        pTxTransducer = xdc_focused_array(stTrans.nNumEle_x, ...
                                          stTrans.nWidth_x, ...
                                          stTrans.nWidth_y, ...
                                          stTrans.nPitch_x-stTrans.nWidth_x, ...
                                          stTrans.nEleFocalDepth, ...
                                          nSubDivNum_x, ...
                                          nSubDivNum_y, ...
                                          aFocalPointPos);
        pRxTransducer = xdc_focused_array(stTrans.nNumEle_x, ...
                                          stTrans.nWidth_x, ...
                                          stTrans.nWidth_y, ...
                                          stTrans.nPitch_x-stTrans.nWidth_x, ...
                                          stTrans.nEleFocalDepth, ...
                                          nSubDivNum_x, ...
                                          nSubDivNum_y, ...
                                          aFocalPointPos);
    elseif strcmp(stTrans.sType,'Convex')
        % Has no elevational focusing lens 
        pTxTransducer = xdc_convex_array (stTrans.nNumEle_x, ...
                                          stTrans.nWidth_x, ...
                                          stTrans.nWidth_y, ...
                                          stTrans.nPitch_x-stTrans.nWidth_x, ...
                                          stTrans.nRadius, ...
                                          nSubDivNum_x, ...
                                          nSubDivNum_y, ...
                                          aFocalPointPos);
        pRxTransducer = xdc_convex_array (stTrans.nNumEle_x, ...
                                          stTrans.nWidth_x, ...
                                          stTrans.nWidth_y, ...
                                          stTrans.nPitch_x-stTrans.nWidth_x, ...
                                          stTrans.nRadius, ...
                                          nSubDivNum_x, ...
                                          nSubDivNum_y, ...
                                          aFocalPointPos);
    else
        error('undefined');
    end
    
    if(bShowFieldIISetUp)
        figure; 
        show_xdc(pTxTransducer); title('Tx Aperture');
    end
    
    %  ??
    xdc_baffle(pTxTransducer,0);      
    xdc_baffle(pRxTransducer,0); 

    %%% Set impulse response and excitation pulse so that the round-trip signal can have a bandwidth of 'nSignalBandwidth'
   
%     %%% Impulse Response (2.2: -6dB BW is 88% when nFc = 5.208 MHz)
%     at2 = 0:dt:2.2*nTc;
%     aTransImpulsResp = sin(2*pi*nFc*at2); 
%     aTransImpulsResp = aTransImpulsResp.*(hanning(max(size(aTransImpulsResp)))'); % Transducer's impulse response
%     xdc_impulse(pTxTransducer,aTransImpulsResp); % Tx impulse response 
%     xdc_impulse(pRxTransducer,aTransImpulsResp); % Rx impulse response 
% 
%     
%     %%% Excitation Pulse
%     at1 = 0:dt:nTxPulseCycle*nTc;
%     nPulseAmp = 1;
%     aPulseSeq  = nPulseAmp*sin(2*pi*nFc*at1);  % excitation pulse
%     xdc_excitation(pTxTransducer,aPulseSeq);  % convol Tx aperture and pulse
% 
%     %%% Impulse Response (2.2: -6dB BW is 88% when nFc = 5.208 MHz)
%     at2 = -nImpulseRespPulseCycle*nTc/2:dt:nImpulseRespPulseCycle*nTc/2;
%     aTransImpulsResp = cos(2*pi*nFc*at2); 
%     aTransImpulsResp = aTransImpulsResp.*(hanning(max(size(aTransImpulsResp)))'); % Transducer's impulse response
%     xdc_impulse(pTxTransducer,aTransImpulsResp); % Tx impulse response 
%     xdc_impulse(pRxTransducer,aTransImpulsResp); % Rx impulse response 

    %%% Impulse Response 
    if isfield(stTrans,'nBW_45dB')
        nBW = stTrans.nBW_45dB/100; % 0~1
        nDecibel = 45; % [dB]
    elseif  isfield(stTrans,'nBW_6dB')
        nBW = stTrans.nBW_6dB/100; % 0~1
        nDecibel = 6; % [dB]
    else
        error('Bandwidth of transducer is undefined: nBW_45dB or nBW_6dB must be clarified in %');
    end
    t_b = 2*sqrt(nDecibel/20*log(10))/(pi*nBW*nFc);% when BW is -nDecibel dB bandwidth
    at2 = (-8/sqrt(2)*t_b):dt:(8/sqrt(2)*t_b); % plot gaussian in a range of +/-8std (std of Guassian = t_b/sqrt(2))
    aEnvelop = exp(-((at2-mean(at2))/t_b).^2);
    aTransImpulsResp = aEnvelop.*sin(2*pi*nFc*at2); 
    xdc_impulse(pTxTransducer,aTransImpulsResp); % Tx impulse response 
    xdc_impulse(pRxTransducer,aTransImpulsResp); % Rx impulse response 
    
    %%% Excitation Pulse
    at1 = -nTxPulseCycle*nTc/2:dt:nTxPulseCycle*nTc/2;
    nPulseAmp = 1;
    aPulseSeq  = nPulseAmp*cos(2*pi*stTx.nFc*at1);  % excitation pulse
    xdc_excitation(pTxTransducer,aPulseSeq);  % convol Tx aperture and pulse

    if(bShowFieldIISetUp)
        % Excitation pulse, Impulse response
        figure;
        subplot(2,3,1)
            plot(at1*1e6,aPulseSeq); title('Exitation Pulse')
        subplot(2,3,2)
            plot(at2*1e6,aTransImpulsResp); title('Transducer`s Impulse response')
        subplot(2,3,3)
            aConvSeq = conv(conv(aPulseSeq,aTransImpulsResp),aTransImpulsResp);
            plot(0:dt*1e6:dt*(length(aConvSeq)-1)*1e6,aConvSeq); title('Expected received signal')
        subplot(2,3,5)
            af = linspace(-1/dt/2,1/dt/2,4097);
            aRTS_f = abs(fftshift(fft(aTransImpulsResp, 4097)));
            aRTS_db = db(aRTS_f/max(aRTS_f)); % in dB scale
            plot(af(2049:end)*1e-6,aRTS_db(2049:end)); xlabel('freq(MHz)'); title('Freq Spectrum of Txdcr Impulse Resp');
            ylim([-150 0]);xlim([0 nFs*1e-6]); grid on; grid minor;
        subplot(2,3,6)
            af = linspace(-1/dt/2,1/dt/2,4097);
            aRTS_f = abs(fftshift(fft(aConvSeq, 4097)));
            aRTS_db = db(aRTS_f/max(aRTS_f)); % in dB scale
            plot(af(2049:end)*1e-6,aRTS_db(2049:end)); xlabel('freq(MHz)'); title('Freq Spectrum of Expected RX signal');
            ylim([-150 0]);xlim([0 nFs*1e-6]); grid on; grid minor;
    end
    nLag = round(((numel(aPulseSeq)+numel(aTransImpulsResp)-1)+numel(aTransImpulsResp)-1)/2); % for received signal
%     nLag = round((numel(aPulseSeq)+numel(aTransImpulsResp)-1)/2); % for pressure field
        
    %%% Take back source information for delay calculation
    nSrcNum = nSubDivNum_x*nSubDivNum_y;
    mSrcPos = GetSrcPos(pTxTransducer, stTrans.nNumEle_x, stTrans.nNumEle_y, nSubDivNum_x, nSubDivNum_y, stTrans.sType, stTrans.nRadius);
    
    %%% num of samples of RF data 
%     nMaxDepth = sqrt(max(stScat.mScatXYZPos(:,3))^2 + (max(stTrans.mElePos(:,1))*2)^2);
%     nSample  = round(nMaxDepth*2/nSoundSpeed*nFs);
    nSample = 2048; disp('!!!!!!RcvData Sample is fixed at 2048!!!!!!!!!');
    vRcvData = zeros(nSample, stTrans.nNumEle_x, stTx.nNum);
    
    for txidx = 1:stTx.nNum
        
        disp(['txidx=' num2str(txidx) '/' num2str(stTx.nNum)]);
        
        %% Tx Aperture & delay 
        % tx angle (CF: scanline angle, PW: plane wave angle)
%         nTxAngle = stTx.aPwAngle_deg(txidx); % [degree]
        
        % aperture
        if isfield(stTx,'mApod')
            aTxAperture = stTx.mApod(txidx,:);
        else
            aTxAperture = TxAptrGen(stTrans, stTx, txidx);
        end
        xdc_apodization(pTxTransducer, 0, aTxAperture);
        
        % delay  
        if strcmp(stTx.sType, 'CF')
            [aTxDelay, nDelayOffset] = TxDelayGen(stTrans, stTx, mSrcPos, nSrcNum, 0, stTx.nNum, txidx, aTxAperture, nSoundSpeed, 'first firing time');
        elseif strcmp(stTx.sType, 'PW')
            [aTxDelay, nDelayOffset] = TxDelayGen(stTrans, stTx, mSrcPos, nSrcNum, stTx.aPwAngle_deg(txidx), 0, 0, aTxAperture, nSoundSpeed, 'first firing time');
        end
        xdc_focus_times(pTxTransducer,0,aTxDelay); % Tx Delay setting
        
        % export
        stTx.aDelayOffset(txidx)      = nDelayOffset; % [sec]
        stTx.mDelay(txidx,:)          = aTxDelay;
        stTx.mAptr(txidx,:)           = aTxAperture;        
        
        %% Rx Delay
        aRxDelay = zeros(1,stTrans.nNumEle_x);
        xdc_center_focus(pRxTransducer,[0 0 0]);% Set the origin for the dynamic focusing line.
        xdc_focus_times(pRxTransducer,0,aRxDelay); % Rx Delay setting
        
        %% Backscattered signal calc %            
        display('Backscattered signal calculation');
        
        if any(aTxAperture)

            %%% Calculate the received signals from a collection of scatterers for all elements
            [mRF, nStartTime] = calc_scat_multi(pTxTransducer, pRxTransducer, stScat.mScatXYZPos, stScat.aScatMag);    

            %%% Zero-padding (so that the first sample is captured at t = 0)
            %%% Modified: replaced padarray with basic MATLAB (no Image Processing Toolbox)
            nPadLen = round(nStartTime/dt);
            % mRF_pad = padarray(mRF,[nPadLen 0],'pre');  % Original (requires Image Processing Toolbox)
            mRF_pad = [zeros(nPadLen, size(mRF,2)); mRF];  % Replacement

            %%% Truncation or Zero-Padding (so that # of samples = RFInfo.nSample)
            nUpSample = nSample*nUpSampFactor;
            nFirstSamIdx = 1 + nLag;
            nLastSamIdx = nFirstSamIdx + nUpSample -1;

            if(size(mRF_pad,1) > nLastSamIdx)
                mRF_trc = mRF_pad(nFirstSamIdx:nLastSamIdx,:);
            elseif(size(mRF_pad,1) < nLastSamIdx)
                %%% Modified: replaced padarray with basic MATLAB (no Image Processing Toolbox)
                % mRF_pad2 = padarray(mRF_pad,[nLastSamIdx-size(mRF_pad,1) 0],'post');  % Original
                mRF_pad2 = [mRF_pad; zeros(nLastSamIdx-size(mRF_pad,1), size(mRF_pad,2))];  % Replacement
                mRF_trc = mRF_pad2(nFirstSamIdx:nLastSamIdx,:);
            end

            %%% Decimate to lower sampling frequency (nFs)
            mRF_deci = mRF_trc(1:nUpSampFactor:end,:);


        else
            mRF_deci = zeros(nSample,nNumRxEle);        
        end

        %%% Save RF data
%         display('Save RF data');
        mRcvData = mRF_deci;


        %%% TGC
%         nTGC_dB = 0.5; % [dB/(MHz*cm)]
%         nSoundSpeed = 1540; % [m/s]
%         aTime = (0:(size(mRcvData,1)-1)).'/nFs + 0;
%         aZ_cm = aTime*nSoundSpeed/2 * 1e2;
%         aTGC = 10.^(nTGC_dB*(nFc/1e6)*aZ_cm/20);
%         mTGC = repmat(aTGC,[1,size(mRcvData,2)]);
%         mRcvData = mRcvData .* mTGC;


        vRcvData(:,:,txidx) = mRcvData;


     
    end
    
    %% Export
    
    stRcvData.sVersion    = sVersion;
    stRcvData.vRcvData    = vRcvData;
%     stRcvData.aDataSize   = size(vRcvData);
    stRcvData.stTx        = stTx;
    stRcvData.nFc         = nFc;
    stRcvData.nFs         = nFs;
    stRcvData.nSoundSpeed = nSoundSpeed;
    stRcvData.stScat      = stScat;
    stRcvData.stTrans     = stTrans;
    stRcvData.sStartAcqTime = sStartAcqTime;


%     save([sDataDir 'RcvData\' 'stRcvData.mat'], 'stRcvData');
%     fid = fopen([sDataDir 'RcvData\' 'vRcvData.bin'],'wb');
%     fwrite(fid, vRcvData, 'float');
%     fclose(fid);

    
end
    



