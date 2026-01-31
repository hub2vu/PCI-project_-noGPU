% Sua Bae
%
% 2019-05-23
%   update at calculating PW LWR region
%        -or+ pitch/2 : to include aperture from end to end (not to reject side scanlines of multibeams)
%
% 2019-05-22
%  RF beamformer for plane wave imaging
%  linear
% 
function RFBeamformer(sRfDir, sBfDir, stRFInfo, stBFpm)

    if ~strcmp(stRFInfo.sVersion,'2.1')
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
    
    % BF information 
    nChannel_x      = stRFInfo.nRcvChannel; % num. of receiving channels
    nSoundSpeed     = stRFInfo.nSoundSpeed; % [m/s]
    nRxFnum         = stBFpm.nRxFnum;
    sRxApodWindow   = stBFpm.sRxApodWindow;  
    nSysOffset_sec = stBFpm.nSysOffset_sec; % to vary the offset [sec] % nSysOffset_sec: residual offset caused by system error (?)
                
    % Demodulation
    stDemod         = stBFpm.stDemod;
    
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
            mDiff = bsxfun(@minus, stBFpm.stG_rt.aT, stTx.mFocalPos(:,5)'); % theta of BF - theta of TX, size: (stG_rt.nTdim x stTx.nNum)
        elseif strcmp(sTransType,'Convex')
            mDiff = bsxfun(@minus, stBFpm.stG.aZ, stTx.mFocalPos(:,3)'); % depth of BF - depth of TX focus, size: (stG.nZdim x stTx.nNum)
        end
        [~, aTxIdx] = min(abs(mDiff),[],2); % data index (txidx) for reconstruction of each BF scanline, size: (stG_rt.nTdim x 1)
    end
    
    mRFBFOut = zeros(size(mP_x));
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
            aScIdx_bf = find(aTxIdx == txidx); % scanline indices for reconstruction using txidx-th data (aScIdx_bf: recon (BF) scanline indices, txidx: tx scanline index)
                                                % For single beam, aScIdx_bf = 1 for txidx = 1
                                                % For 4 multi beam, aScIdx_bf = [1,2,3,4] for txidx = 1
%             display(['    ' num2str(txidx) '-th scanline data (' num2str(aFocalPos(5)) ' deg) is being processed for ' num2str(aThtIdx_bf') '-th BF scanline (' num2str(stBFpm.stG_rt.aT(aThtIdx_bf)') ' deg)']);   
        end
        
        % load data
        load([sRfDir '\mRcvData_tx_' num2str(txidx) '.mat'],'mRcvData');
        if size(mRcvData,2) ~= nChannel_x
            error('mismatched nChannel_x (stRFInfo.nRcvChannel)')
        end
        
        %%%     0) Decimation and Interpolation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isfield(stBFpm, 'nDeciFactor');
            if stBFpm.nDeciFactor >= 1
                mRcvData_deci = mRcvData(1:stBFpm.nDeciFactor:end,:);
                nFs_deci = stRFInfo.nFs/stBFpm.nDeciFactor;
            else
                error('undefined');
            end
        else
            mRcvData_deci = mRcvData;
            nFs_deci = stRFInfo.nFs;
        end
        if isfield(stBFpm, 'nIntpFactor');
            if stBFpm.nIntpFactor >= 1
                mRcvData_intp = Interp_v1(mRcvData_deci, stBFpm.nIntpFactor);
                nFs_intp = nFs_deci*stBFpm.nIntpFactor;
            else                
                error('undefined');
            end
        else
            mRcvData_intp = mRcvData_deci;
            nFs_intp = nFs_deci;
        end
        
        nRFSample = size(mRcvData_intp,1); % num. of samples  
        
% f=figure; 
% nNumFft = 14097;
% subplot(1,2,1); 
%     [af_half, aS_db_half] = myFFT(mRcvData(:,64), nNumFft, stRFInfo.nFs);
%     plot(af_half, aS_db_half); xlim([min(af_half),max(af_half)]); ylim([-50 0]); grid minor; title('original');
% subplot(1,2,2); 
%     [af_half, aS_db_half] = myFFT(mRcvData_intp(:,64), nNumFft, nFs_intp);          
%     % find -6dB width
%     [~, maxidx] = max(aS_db_half);
%     [~, lftidx] = min(abs(aS_db_half(1:maxidx)+6));
%     [~, rgtidx_tmp] = min(abs(aS_db_half(maxidx+1:end)+6));
%     rgtidx = rgtidx_tmp + maxidx;
%     % plot
%     hplot = plot(af_half, aS_db_half); xlim([min(af_half),max(af_half)]); ylim([-50 0]); grid minor; title('decimated&interpolated');
%     % data cursor tip
%     hCursorMode = datacursormode(f);
%     hDatatip(1) = hCursorMode.createDatatip(hplot);
%     hDatatip(1).Position = [af_half(lftidx),aS_db_half(lftidx)];% left
%     hDatatip(2) = hCursorMode.createDatatip(hplot);
%     hDatatip(2).Position = [af_half(rgtidx),aS_db_half(rgtidx)];% right
%     % calc %BW
%     nMaxFreq = FindCenterFreq(mRcvData_intp(:,64), nFs_intp, nNumFft);   
%     disp(['BW = ' num2str( (af_half(rgtidx) - af_half(lftidx) )/nMaxFreq*100 ) '%']);
     
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
            mLogic_PR(:,aScIdx_bf) = 1;
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
        aBfData = zeros(size(aBF_x)); % beamformed data (only within LWR)
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
                aAddress = aTimeDelay*nFs_intp;     
                aAddress = max(min(aAddress, nRFSample-1),0);% saturation, if Add=0, summation X
                aLogic   = ( aAddress > 0 ).* (aAddress < nRFSample-1);

                % 4) spline-interpolation
                bIntp = isfield(stBFpm, 'nIntpFactor');  
                if bIntp % if data is already interpolated, do not interpolate more
                    aRfData = mRcvData_intp(:,cidx);
                    aPicked = aRfData(round(aAddress)+1);
                else
                    aPicked = interpn(0:1:(nRFSample-1), mRcvData_intp(:,cidx), aAddress,'spline');
                end

                % 5) accumulate
                aBfData = aBfData + aPicked.*aLogic.*aApod;  
            end

        end
        
        %%%     3) Compounding  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mRFBFOut(mLogic_PR) = mRFBFOut(mLogic_PR) + aBfData;
        mCompNum(mLogic_PR) = mCompNum(mLogic_PR) + 1;
        if(stBFpm.bSaveEachAngle)
            mRFBFOut_eachAng                = zeros(size(mP_x));
            mRFBFOut_eachAng(mLogic_PR)     = aBfData ; 
            stBfData_b4c.mBfData            = mRFBFOut_eachAng;
            stBfData_b4c.stRFInfo           = stRFInfo;
            stBfData_b4c.stBFpm             = stBFpm;
            save([sBfDir, 'stBfData_b4c_tx_' num2str(txidx)],'stBfData_b4c');
        end
        if(stBFpm.bPlot); axes(axcmp); imagesc(mP_x(1,:),mP_z(:,1),db(abs(mRFBFOut)/max(max(abs(mRFBFOut)))));axis equal; axis tight; caxis([-100 0]);pause(0.1); 
        end;

    
    end
    
    %%%     4) Normalization   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    mRFBFOut = mRFBFOut./mCompNum;
    mRFBFOut(isnan(mRFBFOut)) = 0;
    
    %%   4. Demodulation
    mIQOut = QDM_v2(mRFBFOut, nFs_intp, stRFInfo.nSoundSpeed, stDemod, 1); % without decimation in RFBF
        
    
    %%   4. DSC to catesian grid (only for Convex data)
    switch sTransType
        case 'Linear'
            mBfData = abs(mIQOut);
            % save mBfData
            stBfData.mBfData     = mBfData;
            
        case 'Convex'
            
            % save before DSC
            stBfData.mBfData_rt  = abs(mIQOut); %  (rt grid) 
            
            if stBFpm.bDSC %%%  FROM RT GRID TO XZ GRID
                stG_rt = stBFpm.stG_rt; % BF grid (rt grid)
                stG     = stBFpm.stG; % DSC grid (xz grid)
                mBfData = DSC_rt2xz(abs(mIQOut), stG_rt, stG);
                
                % save after DSC
                stBfData.mBfData     = mBfData;
                
                % plot
                if(stBFpm.bPlot);figure;imagesc(stG.aX,stG.aZ0,db(abs(mBfData)/max(abs(mBfData(:)))));axis equal; axis tight; caxis([-100 0]);title('after dsc'); end;
            end
            
    end
    
    %%   4. Save else
    stBfData.mCompNum    = mCompNum; % when 'Convex', it is on RT grid
    stBfData.stRFInfo    = stRFInfo;
    stBfData.stBFpm      = stBFpm;
    save([sBfDir stBFpm.sBfOutName], 'stBfData');
    display(sprintf('%d/%02d/%d %02d:%02d:%02.2f',clock));

end