% Sua Bae
%
% 2019-07-08
%   update for adding option : stBFpm.bExportIQData
% 2019-05-23
%   update at calculating PW LWR region
%        -or+ pitch/2 : to include aperture from end to end (not to reject side scanlines of multibeams)
% 2019-04-05
%   update for stBFpm.sVersion 2.0->3.0
%    add Interpolation and phase quantization (nQiq(nIntpFactor), nQph)
%
% 2019-03-14
%   update for stBFpm.sVersion 2.0->2.1
%    add stBFpm.nMultiBeam 
%    add stBFpm.sBfOutName
%    add stBFpm.bDSC
%
% 2019-02-07
%   update for 'stRFInfo v2.0' 
%        using 'stRFInfo.nRcvChannel'
%        using 'stBFpm.stDemod' and 'QDM_v2.m'
%   update for convex+CF (E-cube)
%   update for finding leftmost & rightmost element of tx aperture
%        if 'stTx.mAptr' exist, use it! Else make it using TxAptrGen
%   update for stTx.aDelayOffset
%        if stRFInfo.sStartAcqTime = 'first firing time', stTx.aDelayOffset(txidx) must be specified
%   update for 'stRFInfo.aRxAptrPos(txidx)'
%        to compensate rx aperture position when sub-aperture used for RX, ex) 64-channel RX (E-cube)
%   update for dynamic demodulation
%        stBFpm.stDemod: bDynamicDemod
%   update for exclusion of dummy channels of walking aperture from channel-sum
%
% 2017-12-21
%   convex, conventional focusing (simulation)
% 
% 2017-12-19 (LINEAR & CONVEX, PLANE WAVE & DIVERGING WAVE)
%   update for using stG_rt
%   include DSC for convex
%   include diverging wave beamforming
%
% 2017-08-28
%  IQ beamformer for plane wave imaging
%  lienar & convex
% 
function IQBeamformer(sRfDir, sBfDir, stRFInfo, stBFpm)

    if ~(strcmp(stRFInfo.sVersion,'2.1')||strcmp(stRFInfo.sVersion,'2.6'))
        % nRcvChannel should be included
        error('version mismatch');
    end
    if ~strcmp(stBFpm.sVersion,'3.0')
        % nTxOffsetDelay_sec -> nSysOffset_sec
        % add stBFpm.nMultiBeam 
        % add stBFpm.bDSC
        error('version mismatch'); 
    end
    
    %% 1. BF Parameters Download
%     display('    1. BF Parameters Download');
    
    % Transducer
    stTrans             = stRFInfo.stTrans;
    nTransEleNum_x      = stTrans.nNumEle_x; % num of elements in lateral
    nTransElePitch_x    = stTrans.nPitch_x; % [meter]
    nTransRadius        = stTrans.nRadius; % [meter]
    
    % stTx
    stTx            = stRFInfo.stTx;    
                
    % Demodulation
    stDemod         = stBFpm.stDemod;
    
    % IQ beamformer de
    nChannel_x      = stRFInfo.nRcvChannel; % num. of receiving channels
    nSoundSpeed     = stRFInfo.nSoundSpeed; % [m/s]
    
    % BF information 
    nRxFnum         = stBFpm.nRxFnum;
    sRxApodWindow   = stBFpm.sRxApodWindow;  
    nSysOffset_sec = stBFpm.nSysOffset_sec; % to vary the offset [sec] % nSysOffset_sec: residual offset caused by system error (?)

    
    sTransType = stTrans.sType;
    if strcmp(sTransType,'Linear')
        [mP_z, mP_x] = ndgrid(stBFpm.stG.aZ, stBFpm.stG.aX);
    elseif strcmp(sTransType,'Convex')
        [mR, mT] = ndgrid(stBFpm.stG_rt.aR, stBFpm.stG_rt.aT);
        [mP_x, mP_z] = rt2xz(stBFpm.stG_rt, mR, mT);
    else error('undefined');
    end
    

    %%   2. Calculate the Elements Coordinates & Position of origin of scanline 
%     display('    2. Calculate the Elements Coordinates & Position of origin of scanline ');
    % the center of the transducer = the origin (0,0)
    
    %%%    Tx element position     
    switch sTransType
        case 'Linear'
            aEl_x = stTrans.mElePos(:,1)'; % x-pos of elements
            aEl_z = stTrans.mElePos(:,3)'; % z-pos of elements
        case 'Convex'
            aEl_a = stTrans.mElePos(:,4)'; % Transducer element theta position [degree] % size = [nChannel]
            aEl_x = stTrans.mElePos(:,1)';
            aEl_z = stTrans.mElePos(:,3)';   
    end      


    %%   3. Beamforming & compounding
%     display('    3. Beamforming');    
    
    if strcmp(stTx.sType,'CF')
        if strcmp(sTransType,'Linear')
            mDiff = bsxfun(@minus, stBFpm.stG.aX(:), stTx.mFocalPos(:,1)'); % x of BF - x of TX focus, size: (stG.nXdim x stTx.nNum)
        elseif strcmp(sTransType,'Convex')
            mDiff = bsxfun(@minus, stBFpm.stG_rt.aT(:), stTx.mFocalPos(:,5)'); % theta of BF - theta of TX, size: (stG_rt.nTdim x stTx.nNum)
        end
        [~, aTxIdx] = min(abs(mDiff),[],2); % data index (txidx) for reconstruction of each BF scanline, size: (stG_rt.nTdim x 1)
    end
    
    mIQBFOut = zeros(size(mP_x));
    mCompNum = zeros(size(mP_x));
    if(stBFpm.bPlot); fcmp = figure; axcmp = subplot(1,1,1); end;
    for txidx = 1:stTx.nNum
                
        % Tx info 
        if strcmp(stTx.sType,'PW')
            nTxAngle = stTx.aPwAngle_deg(txidx); % [degree]
%             display(['    ' num2str(txidx) '-th angle is being processed: ' num2str(nTxAngle) 'deg']);        
        elseif strcmp(stTx.sType,'DW')
            aVSPos   = stTx.mVSPos(txidx,:); % (x,y,z) of virtual source position  [meter]  
%             display(['    ' num2str(txidx) '-th virtual source data is being processed: x = ' num2str(aVSPos(1)*1e3) ' mm, z = ' num2str(aVSPos(3)*1e3) ' mm']);   
        elseif strcmp(stTx.sType,'CF')
            aFocalPos   = stTx.mFocalPos(txidx,:); % (x,z,z0,r,tht) of virtual source position [meter, meter, meter, meter, deg]
            aThtIdx_bf = find(aTxIdx == txidx); % theta indices for reconstruction using txidx-th data (theta index: recon (BF) scanline index, txidx: tx scanline index)
                                                % For single beam, aThtIdx_bf = 1 for txidx = 1
                                                % For 4 multi beam, aThtIdx_bf = [1,2,3,4] for txidx = 1
%             display(['    ' num2str(txidx) '-th scanline data (' num2str(aFocalPos(5)) ' deg) is being processed for ' num2str(aThtIdx_bf') '-th BF scanline (' num2str(stBFpm.stG_rt.aT(aThtIdx_bf)') ' deg)']);   
        end
        
        % load data
        load([sRfDir '\mRcvData_tx_' num2str(txidx) '.mat']);
        if size(mRcvData,2) ~= nChannel_x
            error('mismatched nChannel_x (stRFInfo.nRcvChannel)')
        end
        
        %%%     0-1) Demodulation and Decimation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mIQData_tmp = QDM_v2(mRcvData, stRFInfo.nFs, stRFInfo.nSoundSpeed, stDemod, stBFpm.nDeciFactor);     

        %%%     0-2) Interpolation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        bIntp = isfield(stBFpm, 'nIntpFactor');        
        if bIntp&&(stBFpm.nIntpFactor>1)
            mIQData = Interp_v1(mIQData_tmp, stBFpm.nIntpFactor);
            nFs_iq = stRFInfo.nFs/stBFpm.nDeciFactor*stBFpm.nIntpFactor;
        else
            mIQData = mIQData_tmp;             
            nFs_iq = stRFInfo.nFs/stBFpm.nDeciFactor;
        end
        nIQSample = size(mIQData,1); % num. of samples  
        
% f=figure; 
% nNumFft = 14097;
% subplot(1,2,1); 
%     [af_half, aS_db_half] = myFFT(mRcvData(:,64), nNumFft, stRFInfo.nFs);
%     plot(af_half, aS_db_half); xlim([min(af_half),max(af_half)]); ylim([-50 0]); grid minor; title('original');
% subplot(1,2,2); 
%     [af, aS] = myFFT(mIQData(:,64), nNumFft, nFs_iq, 'allFreq');          
%     % find -6dB width
%     [~, maxidx] = max(aS);
%     [~, lftidx] = min(abs(aS(1:maxidx)+6));
%     [~, rgtidx_tmp] = min(abs(aS(maxidx+1:end)+6));
%     rgtidx = rgtidx_tmp + maxidx;
%     % plot
%     hplot = plot(af, aS); xlim([min(af),max(af)]); ylim([-50 0]); grid minor; title('decimated&interpolated');
%     % data cursor tip
%     hCursorMode = datacursormode(f);
%     hDatatip(1) = hCursorMode.createDatatip(hplot);
%     hDatatip(1).Position = [af(lftidx),aS(lftidx)];% left
%     hDatatip(2) = hCursorMode.createDatatip(hplot);
%     hDatatip(2).Position = [af(rgtidx),aS(rgtidx)];% right
%     % calc %BW
%     disp(['BW = ' num2str( (af(rgtidx) - af(lftidx) )/stRFInfo.nFc*100 ) '%']);

        %%%     1) Select BF points reconstructed by using the tidx-th data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % leftmost and rightmost active elements
        if isfield(stTx,'mAptr')
            if size(stTx.mAptr) == [stTx.nNum, stTrans.nNumEle_x]
                aAptr = stTx.mAptr(txidx,:);
            else error('stTx.mAptr: dimension mismatch'); 
            end
        else
            aAptr = TxAptrGen(stTrans, stTx, txidx);
        end
        lftmidx = find(aAptr,1,'first'); 
        rftmidx = find(aAptr,1,'last');
        
        % find beamforming point
        if strcmp(stTx.sType, 'CF')
            % Select Beamforming Points
            mLogic_PR = zeros(size(mP_x,1),size(mP_x,2));
            mLogic_PR(:,aThtIdx_bf) = 1;
            mLogic_PR = logical(mLogic_PR);
            
        elseif strcmp(stTx.sType, 'PW')||strcmp(stTx.sType, 'DW')      
            nLftm_x = aEl_x(lftmidx)-stTrans.nPitch_x/2; nRtm_x  = aEl_x(rftmidx)+stTrans.nPitch_x/2; % -or+ pitch/2 : to include aperture from end to end (not to reject side scanlines of multibeams)
            nLftm_z = aEl_z(lftmidx); nRtm_z  = aEl_z(rftmidx);       
            if strcmp(stTx.sType,'PW')
                % plane wave propagation region
                mLogic_PR = (mP_x >= (nLftm_x + (mP_z-nLftm_z)*tand(nTxAngle))) & (mP_x <= (nRtm_x + (mP_z-nRtm_z)*tand(nTxAngle)));
            elseif strcmp(stTx.sType,'DW')
                nAng_L  = atand( (nLftm_x-aVSPos(1))./(nLftm_z-aVSPos(3)) );% left diverging angle
                nAng_R  = atand( (nRtm_x-aVSPos(1))./(nRtm_z-aVSPos(3)) );% left diverging angle
                % diverging wave propagation region
                mAng_vs2bf = atand( (mP_x-aVSPos(1))./(mP_z-aVSPos(3)) ); % angle between virtual source and BF point   
                mLogic_PR = (mAng_vs2bf >= nAng_L) & (mAng_vs2bf <= nAng_R);
            end

        else
            error('undefined');
        end

        % Select Beamforming Points
        aBF_x = mP_x(mLogic_PR);
        aBF_z = mP_z(mLogic_PR);
        
        %%%     2) Beamforming  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Tx delay for beamforming
        if strcmp(stTx.sType,'PW')
            aTxDist = aBF_x*sind(nTxAngle) + aBF_z*cosd(nTxAngle); % [meter] Distance from the PW at origin to target [meter] 
        elseif strcmp(stTx.sType,'DW')
            aTxDist = sqrt((aVSPos(1)-aBF_x).*(aVSPos(1)-aBF_x) + (aVSPos(3)-aBF_z).*(aVSPos(3)-aBF_z)); % [meter] Distance from the virtual source to target [meter] 
        elseif strcmp(stTx.sType,'CF')
            % multibeam delay type: CASE1 (the same tx delay for all multibeams)
            switch sTransType
                case 'Linear'
                    aTxDist = aBF_z; % [meter] Distance from the x-axis to target [meter] 
                case 'Convex'
                    aTxDist = sqrt(aBF_x.^2 + aBF_z.^2); % [meter] Distance from the origin to target [meter] 
            end
        end
        aTxDelay = aTxDist/nSoundSpeed; % [sec]
        if strcmp(stRFInfo.sStartAcqTime,'placed at origin')
            aTxDelay = aTxDelay; % do  nothing
        elseif strcmp(stRFInfo.sStartAcqTime,'first firing time')
            % offset (time difference) between the time when the wave was virtually at the origin and the time when the first element was fired
            aTxDelay = aTxDelay - stTx.aDelayOffset(txidx); % [sec]
        else
            error('undefined');
        end
        
        % Rx Aperture size
        switch sTransType
            case 'Linear'
                aAptSize_m = aBF_z/nRxFnum + nTransElePitch_x; % plus nTransElePitch_x to prevent nApod being NaN
            case 'Convex'
                aRadial_distance_from_Txdcr = sqrt(aBF_x.^2 + aBF_z.^2) - nTransRadius; % aRadial_distance_from_Txdcr to each BF point
                aAptSize_m = aRadial_distance_from_Txdcr/nRxFnum + nTransElePitch_x; % plus nTransElePitch_x to prevent nAperture being NaN
                % nMaxAptSize_m = (30/180*pi)*nTransRadius*2; % maximum aperture size is limited by acceptance angle
                nMaxAptSize_m = nTransRadius*sind(15)*2; % maximum aperture size is limited by acceptance angle
                aAptSize_m = min(aAptSize_m, nMaxAptSize_m);
        end        
        
        % Rx channel position
        aChanIdx_x = -(nChannel_x-1)/2 : 1 : (nChannel_x-1)/2;
        if nChannel_x ~= nTransEleNum_x % if not full aperture
            aChanIdx_x = aChanIdx_x + stRFInfo.stTx.aRxAptrPos(txidx); % compensate rx aperture position
        end
        switch sTransType
            case 'Linear'
                aCh_x = aChanIdx_x*nTransElePitch_x; % x-pos of elements
                aCh_z = zeros(size(aCh_x)); % z-pos of elements
            case 'Convex'
                aCh_a = aChanIdx_x*nTransElePitch_x/nTransRadius/pi*180; % Transducer element theta position [degree] % size = [nChannel]
                aCh_x = nTransRadius*sind(aCh_a);
                aCh_z = nTransRadius*cosd(aCh_a);   
        end

        % Interpolation & summation
        aInph = zeros(size(aBF_x)); % beamformed data (only within LWR)
        aQuad = zeros(size(aBF_x)); % beamformed data (only within LWR)
        for cidx = 1 : nChannel_x
            
            if abs(aChanIdx_x(cidx)) <= (nTransEleNum_x-1)/2 % to exclude dummy channels of walking aperture

                % 1) Apodization window
                if(nRxFnum~=0)
                    switch sTransType
                        case 'Linear'
                            aDistance_m = abs(aBF_x - aCh_x(cidx));  % from x-pos of BF point to each element      
                        case 'Phased'
                            aDistance_m = abs(aCh_x(cidx)*ones(numel(aBF_data),1));  % from center to each element       
                        case 'Convex'
                            aBF_a = atand(aBF_x./aBF_z);
                            aDistance_m = abs(nTransRadius*sind(aCh_a(cidx) - aBF_a));  % from specific point (intersection of scanline and array) to each element              
                    end

                    aApod = ApodGen(aDistance_m,aAptSize_m,sRxApodWindow);
    %                 idx_a_certain_BFPoint = 500;
    %                 aApod_for_a_certain_BFpoint(cidx) = aApod(idx_a_certain_BFPoint);
                else
                    aApod = ones(size(aAptSize_m));
                end            

                % 2) Rx delay, Round trip delay
                aRxDelay = sqrt( (aBF_x - aCh_x(cidx)).^2 + (aBF_z - aCh_z(cidx)).^2 ) /nSoundSpeed; % [sec]
                aTimeDelay = aTxDelay + aRxDelay + nSysOffset_sec; % [sec]  % nSysOffset_sec: residual offset caused by system error (?)



                % 3) calc address
                aAddress = aTimeDelay*nFs_iq;     
                aAddress = max(min(aAddress, nIQSample-1),0);% saturation, if Add=0, summation X
                aLogic   = ( aAddress > 0 ).* (aAddress < nIQSample-1);

                % 4) spline-interpolation
                bIntp = isfield(stBFpm, 'nIntpFactor');  
                if bIntp % if data is already interpolated, do not interpolate more
                    aInphChData = real(mIQData(:,cidx));
                    aQuadChData = imag(mIQData(:,cidx));
                    aInph_intp = aInphChData(round(aAddress)+1);
                    aQuad_intp = aQuadChData(round(aAddress)+1);
                else
                    aInph_intp = interpn(0:1:(nIQSample-1), real(mIQData(:,cidx)), aAddress,'spline');
                    aQuad_intp = interpn(0:1:(nIQSample-1), imag(mIQData(:,cidx)), aAddress,'spline');
                end

                % 5) compansate phase (phase = 2*pi*f0*delay)
                if stDemod.bDynamic
                    nFs = stRFInfo.nFs;
                    aK = aTimeDelay*nFs; % sample index before decimation (because demodulation was applied before decimation)
                    aCompen_phase = pi*stDemod.nDemodFreq_slop*nSoundSpeed/nFs/nFs*aK.*(aK-1)/2 + aK*2*pi*stDemod.nDemodFreq_intercept/nFs;
                else
                    % aCompen_phase = (aTimeDelay-2*aBF_z/nSoundSpeed)*2*pi*stDemod.nDemodFreq; %[sec]->[rad]
                    aCompen_phase   = (aTimeDelay)*2*pi*stDemod.nDemodFreq; %[sec]->[rad]
                end
                if isfield(stBFpm,'nQph') % if nQph exists, quantize the phase!
                    nStep = 2*pi/stBFpm.nQph; % quantization step [rad]
                    aCompen_phase = round(aCompen_phase/nStep)*nStep; % quantized phase [rad]
                end
                aInph_cps = aInph_intp.*cos(aCompen_phase) - aQuad_intp.*sin(aCompen_phase);
                aQuad_cps = aInph_intp.*sin(aCompen_phase) + aQuad_intp.*cos(aCompen_phase);

                % 6) accumulate
                aInph = aInph + aInph_cps.*aLogic.*aApod;  
                aQuad = aQuad + aQuad_cps.*aLogic.*aApod;   
            end

        end
        
        %%%     3) Compounding  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mIQBFOut(mLogic_PR) = mIQBFOut(mLogic_PR) + (aInph + 1j*aQuad);
        mCompNum(mLogic_PR) = mCompNum(mLogic_PR) + 1;
        if(stBFpm.bSaveEachAngle)
            mIQBFOut_eachAng                = zeros(size(mP_x));
            mIQBFOut_eachAng(mLogic_PR)    = (aInph + 1j*aQuad); 
            stBfData_b4c.mBfData            = mIQBFOut_eachAng;
            stBfData_b4c.stRFInfo            = stRFInfo;
            stBfData_b4c.stBFpm              = stBFpm;
            save([sBfDir, 'stBfData_b4c_tx_' num2str(txidx)],'stBfData_b4c');
        end
        if(stBFpm.bPlot);  axes(axcmp); imagesc(mP_x(1,:),mP_z(:,1),db(abs(mIQBFOut)/max(max(abs(mIQBFOut)))));axis equal; axis tight; caxis([-100 0]);pause(0.1); 
        end;

    
    end
    
    %%%     4) Normalization   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    mIQBFOut = mIQBFOut;%./mCompNum;
    mIQBFOut(isnan(mIQBFOut)) = 0;
    
    
    %%   4. DSC to catesian grid (only for Convex data)
    switch sTransType
        case 'Linear'
            mBfData = abs(mIQBFOut);
            % save mBfData
            stBfData.mBfData     = mBfData;
            
        case 'Convex'
            
            % save before DSC
            stBfData.mBfData_rt  = abs(mIQBFOut); %  (rt grid) 
            
            if stBFpm.bDSC %%%  FROM RT GRID TO XZ GRID
                stG_rt = stBFpm.stG_rt; % BF grid (rt grid)
                stG     = stBFpm.stG; % DSC grid (xz grid)
                mBfData = DSC_rt2xz(abs(mIQBFOut), stG_rt, stG);
                
                % save after DSC
                stBfData.mBfData     = mBfData;
                
                % plot
                if(stBFpm.bPlot);figure;imagesc(stG.aX,stG.aZ0,db(abs(mBfData)/max(abs(mBfData(:)))));axis equal; axis tight; caxis([-100 0]);title('after dsc'); end;
            end
            
    end
    
    %%   4. Save else
    if isfield(stBFpm,'bExportIQData')
        if stBFpm.bExportIQData
            stBfData.mIQData = mIQData;
            stBfData.mIQBFOut = mIQBFOut;
        end
    end
    stBfData.mCompNum    = mCompNum; % when 'Convex', it is on RT grid
    stBfData.stRFInfo    = stRFInfo;
    stBfData.stBFpm      = stBFpm;
    save([sBfDir stBFpm.sBfOutName], 'stBfData');
    display(sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock));

end