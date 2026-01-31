function mBfData = IQBeamformer_GPU(vRcvData, stRFInfo, stBFpm, bUseGPU)
% IQBeamformer_GPU - GPU 가속 IQ 빔포머
%
% 입력:
%   vRcvData  - RF 데이터 [nSample x nChannel x nTx]
%   stRFInfo  - RF 정보 구조체
%   stBFpm    - 빔포밍 파라미터 구조체
%   bUseGPU   - GPU 사용 여부 (true/false)
%
% 출력:
%   mBfData   - 빔포밍된 이미지 [nZ x nX]
%
% 2025-01-31 - GPU 최적화 버전

    %% 파라미터 추출
    stTrans = stRFInfo.stTrans;
    stTx = stRFInfo.stTx;
    stDemod = stBFpm.stDemod;
    
    nChannel_x = stRFInfo.nRcvChannel;
    nTransEleNum_x = stTrans.nNumEle_x;
    nTransElePitch_x = stTrans.nPitch_x;
    nSoundSpeed = stRFInfo.nSoundSpeed;
    nRxFnum = stBFpm.nRxFnum;
    sRxApodWindow = stBFpm.sRxApodWindow;
    nSysOffset_sec = stBFpm.nSysOffset_sec;
    
    %% 빔포밍 그리드 생성
    [mP_z, mP_x] = ndgrid(stBFpm.stG.aZ, stBFpm.stG.aX);
    
    %% 채널 위치 계산
    aChanIdx_x = -(nChannel_x-1)/2 : 1 : (nChannel_x-1)/2;
    aCh_x = aChanIdx_x * nTransElePitch_x;
    aCh_z = zeros(size(aCh_x));
    
    %% Element 위치
    aEl_x = stTrans.mElePos(:,1)';
    aEl_z = stTrans.mElePos(:,3)';
    
    %% GPU로 데이터 전송
    if bUseGPU
        mP_x = gpuArray(single(mP_x));
        mP_z = gpuArray(single(mP_z));
        aCh_x = gpuArray(single(aCh_x));
        aCh_z = gpuArray(single(aCh_z));
    end
    
    %% 빔포밍 출력 초기화
    if bUseGPU
        mIQBFOut = gpuArray(complex(zeros(size(mP_x), 'single')));
    else
        mIQBFOut = complex(zeros(size(mP_x)));
    end
    
    %% 각 Tx angle에 대해 빔포밍
    for txidx = 1:stTx.nNum
        
        % RF 데이터 로드
        mRcvData = vRcvData(:,:,txidx);
        
        %% IQ Demodulation
        mIQData_tmp = QDM_v2(mRcvData, stRFInfo.nFs, stRFInfo.nSoundSpeed, stDemod, stBFpm.nDeciFactor);
        
        %% Interpolation
        if isfield(stBFpm, 'nIntpFactor') && (stBFpm.nIntpFactor > 1)
            mIQData = Interp_v1(mIQData_tmp, stBFpm.nIntpFactor);
            nFs_iq = stRFInfo.nFs / stBFpm.nDeciFactor * stBFpm.nIntpFactor;
        else
            mIQData = mIQData_tmp;
            nFs_iq = stRFInfo.nFs / stBFpm.nDeciFactor;
        end
        nIQSample = size(mIQData, 1);
        
        %% GPU로 IQ 데이터 전송
        if bUseGPU
            mIQData = gpuArray(single(mIQData));
        end
        
        %% Tx aperture 정보
        if isfield(stTx, 'mAptr')
            aAptr = stTx.mAptr(txidx,:);
        else
            aAptr = ones(1, nTransEleNum_x);
        end
        lftmidx = find(aAptr, 1, 'first');
        rftmidx = find(aAptr, 1, 'last');
        
        %% PW 전파 영역 계산
        nTxAngle = stTx.aPwAngle_deg(txidx);
        nLftm_x = aEl_x(lftmidx) - nTransElePitch_x/2;
        nRtm_x = aEl_x(rftmidx) + nTransElePitch_x/2;
        nLftm_z = aEl_z(lftmidx);
        nRtm_z = aEl_z(rftmidx);
        
        mLogic_PR = (mP_x >= (nLftm_x + (mP_z-nLftm_z)*tand(nTxAngle))) & ...
                    (mP_x <= (nRtm_x + (mP_z-nRtm_z)*tand(nTxAngle)));
        
        %% Tx delay 계산
        aTxDelay = mP_z / nSoundSpeed + mP_x * sind(nTxAngle) / nSoundSpeed;
        
        if strcmp(stRFInfo.sStartAcqTime, 'first firing time')
            aTxDelay = aTxDelay - stTx.aDelayOffset(txidx);
        end
        
        %% Rx aperture size
        aBF_z = mP_z(mLogic_PR);
        aBF_x = mP_x(mLogic_PR);
        aAptSize_m = aBF_z / nRxFnum + nTransElePitch_x;
        
        %% GPU 가속 채널 합산 (벡터화)
        if bUseGPU
            [aInph, aQuad] = beamform_channels_gpu(mIQData, aBF_x, aBF_z, aCh_x, aCh_z, ...
                aTxDelay(mLogic_PR), aAptSize_m, nFs_iq, nIQSample, nSoundSpeed, ...
                stDemod.nDemodFreq, sRxApodWindow, nSysOffset_sec, aChanIdx_x, nTransEleNum_x);
        else
            [aInph, aQuad] = beamform_channels_cpu(mIQData, aBF_x, aBF_z, aCh_x, aCh_z, ...
                aTxDelay(mLogic_PR), aAptSize_m, nFs_iq, nIQSample, nSoundSpeed, ...
                stDemod.nDemodFreq, sRxApodWindow, nSysOffset_sec, aChanIdx_x, nTransEleNum_x);
        end
        
        %% Compounding
        mIQBFOut(mLogic_PR) = mIQBFOut(mLogic_PR) + complex(aInph, aQuad);
        
    end
    
    %% 결과 처리
    mIQBFOut(isnan(mIQBFOut)) = 0;
    mBfData = abs(mIQBFOut);
    
    %% GPU에서 CPU로 전송
    if bUseGPU
        mBfData = gather(mBfData);
    end
    
end

%% ========== GPU 채널 합산 함수 ==========
function [aInph, aQuad] = beamform_channels_gpu(mIQData, aBF_x, aBF_z, aCh_x, aCh_z, ...
    aTxDelay, aAptSize_m, nFs_iq, nIQSample, nSoundSpeed, nDemodFreq, sRxApodWindow, ...
    nSysOffset_sec, aChanIdx_x, nTransEleNum_x)

    nChannel = length(aCh_x);
    nBFPoints = length(aBF_x);
    
    % 모든 채널에 대해 동시에 계산 (벡터화)
    % [nBFPoints x nChannel] 크기의 행렬로 확장
    mBF_x = repmat(aBF_x(:), 1, nChannel);
    mBF_z = repmat(aBF_z(:), 1, nChannel);
    mCh_x = repmat(aCh_x(:)', nBFPoints, 1);
    mCh_z = repmat(aCh_z(:)', nBFPoints, 1);
    mTxDelay = repmat(aTxDelay(:), 1, nChannel);
    mAptSize = repmat(aAptSize_m(:), 1, nChannel);
    
    % Rx delay 계산 (모든 채널 동시)
    mRxDelay = sqrt((mBF_x - mCh_x).^2 + (mBF_z - mCh_z).^2) / nSoundSpeed;
    mTimeDelay = mTxDelay + mRxDelay + nSysOffset_sec;
    
    % Address 계산
    mAddress = mTimeDelay * nFs_iq;
    mAddress = max(min(mAddress, nIQSample-1), 0);
    mLogic = (mAddress > 0) .* (mAddress < nIQSample-1);
    
    % Apodization
    mDistance = abs(mBF_x - mCh_x);
    mApod = ApodGen_vectorized(mDistance, mAptSize, sRxApodWindow);
    
    % 유효 채널 마스크
    mChanMask = abs(repmat(aChanIdx_x(:)', nBFPoints, 1)) <= (nTransEleNum_x-1)/2;
    
    % IQ 데이터 인터폴레이션 (nearest neighbor for speed)
    mAddressIdx = round(mAddress) + 1;
    mAddressIdx = max(min(mAddressIdx, nIQSample), 1);
    
    % 각 채널의 IQ 데이터 추출
    mInphData = real(mIQData);
    mQuadData = imag(mIQData);
    
    % 인덱싱을 위한 선형 인덱스 계산
    [nRows, nCols] = size(mAddressIdx);
    mColIdx = repmat(1:nCols, nRows, 1);
    mLinIdx = sub2ind(size(mIQData), mAddressIdx, mColIdx);
    
    mInph_intp = mInphData(mLinIdx);
    mQuad_intp = mQuadData(mLinIdx);
    
    % Phase compensation
    mCompen_phase = mTimeDelay * 2 * pi * nDemodFreq;
    mInph_cps = mInph_intp .* cos(mCompen_phase) - mQuad_intp .* sin(mCompen_phase);
    mQuad_cps = mInph_intp .* sin(mCompen_phase) + mQuad_intp .* cos(mCompen_phase);
    
    % 채널 합산
    mWeight = mLogic .* mApod .* mChanMask;
    aInph = sum(mInph_cps .* mWeight, 2);
    aQuad = sum(mQuad_cps .* mWeight, 2);
    
end

%% ========== CPU 채널 합산 함수 ==========
function [aInph, aQuad] = beamform_channels_cpu(mIQData, aBF_x, aBF_z, aCh_x, aCh_z, ...
    aTxDelay, aAptSize_m, nFs_iq, nIQSample, nSoundSpeed, nDemodFreq, sRxApodWindow, ...
    nSysOffset_sec, aChanIdx_x, nTransEleNum_x)

    nChannel = length(aCh_x);
    aInph = zeros(size(aBF_x));
    aQuad = zeros(size(aBF_x));
    
    for cidx = 1:nChannel
        if abs(aChanIdx_x(cidx)) <= (nTransEleNum_x-1)/2
            
            % Apodization
            aDistance_m = abs(aBF_x - aCh_x(cidx));
            aApod = ApodGen(aDistance_m, aAptSize_m, sRxApodWindow);
            
            % Rx delay
            aRxDelay = sqrt((aBF_x - aCh_x(cidx)).^2 + (aBF_z - aCh_z(cidx)).^2) / nSoundSpeed;
            aTimeDelay = aTxDelay + aRxDelay + nSysOffset_sec;
            
            % Address
            aAddress = aTimeDelay * nFs_iq;
            aAddress = max(min(aAddress, nIQSample-1), 0);
            aLogic = (aAddress > 0) .* (aAddress < nIQSample-1);
            
            % Interpolation
            aInphChData = real(mIQData(:, cidx));
            aQuadChData = imag(mIQData(:, cidx));
            aInph_intp = aInphChData(round(aAddress)+1);
            aQuad_intp = aQuadChData(round(aAddress)+1);
            
            % Phase compensation
            aCompen_phase = aTimeDelay * 2 * pi * nDemodFreq;
            aInph_cps = aInph_intp .* cos(aCompen_phase) - aQuad_intp .* sin(aCompen_phase);
            aQuad_cps = aInph_intp .* sin(aCompen_phase) + aQuad_intp .* cos(aCompen_phase);
            
            % Accumulate
            aInph = aInph + aInph_cps .* aLogic .* aApod;
            aQuad = aQuad + aQuad_cps .* aLogic .* aApod;
            
        end
    end
    
end

%% ========== 벡터화된 Apodization 함수 ==========
function mApod = ApodGen_vectorized(mDistance, mAptSize, sWindow)
    
    mApod = zeros(size(mDistance));
    mValid = mDistance <= mAptSize/2;
    
    switch sWindow
        case 'boxcar'
            mApod(mValid) = 1;
        case 'hanning'
            mApod(mValid) = 0.5 * (1 + cos(2*pi*mDistance(mValid)./mAptSize(mValid)));
        case 'hamming'
            mApod(mValid) = 0.54 + 0.46 * cos(2*pi*mDistance(mValid)./mAptSize(mValid));
        case 'tukey50'
            mRatio = mDistance ./ mAptSize;
            mTukey = mRatio <= 0.25;
            mApod(mValid & mTukey) = 1;
            mApod(mValid & ~mTukey) = 0.5 * (1 + cos(2*pi*(mRatio(mValid & ~mTukey) - 0.25)/0.5));
        otherwise
            mApod(mValid) = 1;
    end
    
end
