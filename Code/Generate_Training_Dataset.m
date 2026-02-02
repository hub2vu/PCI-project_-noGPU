% Generate_Training_Dataset.m
% 
% U-Net 학습용 데이터셋 생성 코드
% - 128ch (Ground Truth) / 64ch (Input) 동시 생성
% - 같은 시뮬레이션 조건에서 생성되어 매칭됨
% - GPU 가속 빔포밍 적용
% - 랜덤 위치/개수의 점 타겟
%
% 2025-01-31
%
clear; close all;
addpath(genpath('src'));

%% ========== 설정 ==========
nDatasetSize = 100;             % 생성할 데이터 개수
sDatasetDir = '..\dataset\';    % 저장 경로

% 시뮬레이션 파라미터
nMinScatterers = 1;             % 최소 산란체 개수
nMaxScatterers = 10;            % 최대 산란체 개수
nMinDepth_mm = 10;              % 최소 깊이 (mm)
nMaxDepth_mm = 60;              % 최대 깊이 (mm)
nMinLateral_mm = -15;           % 최소 lateral 위치 (mm)
nMaxLateral_mm = 15;            % 최대 lateral 위치 (mm)

% PW 파라미터
ang_deg = 4;                    % max PW angle (deg)
nTxNum = 21;                    % number of PW angles

% 빔포밍 파라미터
nMultiBeam = 1;
nQiq = 2;
nQph = 0;

% GPU 사용 여부 확인
bUseGPU = false;
if gpuDeviceCount > 0
    try
        % CUDA 이후 버전 호환성 활성화 시도
        parallel.gpu.enableCUDAForwardCompatibility(true);
        gpuDevice(1);
        bUseGPU = true;
        disp('GPU detected. Using GPU acceleration (forward compatibility mode).');
    catch ME
        % GPU 사용 실패 시 CPU 모드로 전환
        warning('GPU initialization failed: %s', ME.message);
        disp('Falling back to CPU mode.');
        bUseGPU = false;
    end
else
    disp('No GPU detected. Using CPU.');
end

%% ========== 폴더 생성 ==========
if ~exist(sDatasetDir, 'dir')
    mkdir(sDatasetDir);
end
if ~exist([sDatasetDir 'bf_128ch\'], 'dir')
    mkdir([sDatasetDir 'bf_128ch\']);
end
if ~exist([sDatasetDir 'bf_64ch\'], 'dir')
    mkdir([sDatasetDir 'bf_64ch\']);
end
if ~exist([sDatasetDir 'metadata\'], 'dir')
    mkdir([sDatasetDir 'metadata\']);
end

%% ========== Transducer 로드 ==========
load('stTrans_L7-4.mat');
stTrans.nBW_6dB = 100; % [%] -6dB bandwidth

% 64ch용 트랜스듀서 정보 생성
aChIdx_64 = 1:2:128;  % 홀수 채널 선택
stTrans_64ch = stTrans;
if isfield(stTrans, 'mElePos')
    stTrans_64ch.mElePos = stTrans.mElePos(aChIdx_64, :);
end
if isfield(stTrans, 'nNumEle')
    stTrans_64ch.nNumEle = 64;
end
if isfield(stTrans, 'nNumEle_x')
    stTrans_64ch.nNumEle_x = 64;
end
if isfield(stTrans, 'nPitch_x')
    stTrans_64ch.nPitch_x = stTrans.nPitch_x * 2;
end

%% ========== 데이터 생성 루프 ==========
fprintf('Generating %d training samples...\n', nDatasetSize);
tStart = tic;

for dataIdx = 1:nDatasetSize
    
    fprintf('\n[%d/%d] Generating sample...\n', dataIdx, nDatasetSize);
    tSample = tic;
    
    %% 1. 랜덤 산란체 생성
    nNumScat = randi([nMinScatterers, nMaxScatterers]);
    
    % 랜덤 위치 생성 (x, y, z)
    aX_mm = nMinLateral_mm + (nMaxLateral_mm - nMinLateral_mm) * rand(nNumScat, 1);
    aY_mm = zeros(nNumScat, 1);  % elevation은 0으로 고정
    aZ_mm = nMinDepth_mm + (nMaxDepth_mm - nMinDepth_mm) * rand(nNumScat, 1);
    
    stScat.mScatXYZPos = [aX_mm*1e-3, aY_mm*1e-3, aZ_mm*1e-3];
    stScat.aScatMag = 1e26 * ones(nNumScat, 1);
    
    %% 2. Tx 설정
    stTx.sType = 'PW';
    stTx.nNum = nTxNum;
    stTx.aPwAngle_deg = linspace(-ang_deg, ang_deg, stTx.nNum);
    stTx.nFc = stTrans.nFc;
    stTx.nTxPulseCycle = 1;
    
    %% 3. RF 데이터 생성 (128ch)
    fprintf('  Generating RF data (128ch)...\n');
    bPlot = 0;  % 플롯 비활성화
    stRcvData = GenRFdata(stTrans, stTx, stScat, bPlot);
    
    %% 4. 64ch 데이터 추출
    vRcvData_128ch = stRcvData.vRcvData;
    vRcvData_64ch = vRcvData_128ch(:, aChIdx_64, :);
    
    %% 5. stRFInfo 생성 (128ch, 64ch)
    % 128ch
    stRFInfo_128ch.sVersion = stRcvData.sVersion;
    stRFInfo_128ch.stTrans = stTrans;
    stRFInfo_128ch.nRcvChannel = 128;
    stRFInfo_128ch.nSoundSpeed = stRcvData.nSoundSpeed;
    stRFInfo_128ch.nSample = size(vRcvData_128ch, 1);
    stRFInfo_128ch.nFs = stRcvData.nFs;
    stRFInfo_128ch.nFc = stRcvData.nFc;
    stRFInfo_128ch.stTx = stRcvData.stTx;
    stRFInfo_128ch.stScat = stRcvData.stScat;
    stRFInfo_128ch.sStartAcqTime = stRcvData.sStartAcqTime;
    
    % 64ch
    stRFInfo_64ch = stRFInfo_128ch;
    stRFInfo_64ch.stTrans = stTrans_64ch;
    stRFInfo_64ch.nRcvChannel = 64;
    % stTx 정보도 64ch에 맞게 수정
    stTx_64ch = stRcvData.stTx;
    if isfield(stRcvData.stTx, 'mAptr')
        stTx_64ch.mAptr = stRcvData.stTx.mAptr(:, aChIdx_64);
    end
    if isfield(stRcvData.stTx, 'mDelay')
        stTx_64ch.mDelay = stRcvData.stTx.mDelay(:, aChIdx_64);
    end
    stTx_64ch.aRxAptrPos = zeros(1, stRcvData.stTx.nNum);
    stRFInfo_64ch.stTx = stTx_64ch;
    
    %% 6. GPU 가속 빔포밍
    fprintf('  Beamforming 128ch (GPU)...\n');
    stBFpm_128ch = SetBFParam_simul('IQBF', stRFInfo_128ch, nMultiBeam, nQiq, nQph);
    stBFpm_128ch.bPlot = 0;
    mBfData_128ch = IQBeamformer_GPU(vRcvData_128ch, stRFInfo_128ch, stBFpm_128ch, bUseGPU);
    
    fprintf('  Beamforming 64ch (GPU)...\n');
    stBFpm_64ch = SetBFParam_simul('IQBF', stRFInfo_64ch, nMultiBeam, nQiq, nQph);
    stBFpm_64ch.bPlot = 0;
    mBfData_64ch = IQBeamformer_GPU(vRcvData_64ch, stRFInfo_64ch, stBFpm_64ch, bUseGPU);
    
    %% 7. 데이터 정규화 및 저장
    % 정규화 (0~1 범위)
    mBfData_128ch_norm = abs(mBfData_128ch) / max(abs(mBfData_128ch(:)));
    mBfData_64ch_norm = abs(mBfData_64ch) / max(abs(mBfData_64ch(:)));
    
    % 파일명 생성 (매칭을 위한 동일 ID)
    sDataID = sprintf('%06d', dataIdx);
    
    % BF 데이터 저장
    sFile_128ch = [sDatasetDir 'bf_128ch\bf_' sDataID '.mat'];
    sFile_64ch = [sDatasetDir 'bf_64ch\bf_' sDataID '.mat'];
    
    save(sFile_128ch, 'mBfData_128ch_norm', '-v7.3');
    save(sFile_64ch, 'mBfData_64ch_norm', '-v7.3');
    
    % 메타데이터 저장 (산란체 정보, 시뮬레이션 조건)
    stMetadata.dataID = sDataID;
    stMetadata.nNumScatterers = nNumScat;
    stMetadata.mScatPos_mm = [aX_mm, aY_mm, aZ_mm];
    stMetadata.stTx = stTx;
    stMetadata.ang_deg = ang_deg;
    stMetadata.stBFpm_128ch = stBFpm_128ch;
    stMetadata.stBFpm_64ch = stBFpm_64ch;
    stMetadata.aX_128ch = stBFpm_128ch.stG.aX;
    stMetadata.aZ_128ch = stBFpm_128ch.stG.aZ;
    stMetadata.aX_64ch = stBFpm_64ch.stG.aX;
    stMetadata.aZ_64ch = stBFpm_64ch.stG.aZ;
    stMetadata.timestamp = datetime('now');
    
    sFile_meta = [sDatasetDir 'metadata\meta_' sDataID '.mat'];
    save(sFile_meta, 'stMetadata');
    
    tElapsed = toc(tSample);
    fprintf('  Sample %d completed in %.2f sec\n', dataIdx, tElapsed);
    
end

%% ========== 완료 ==========
tTotal = toc(tStart);
fprintf('\n========================================\n');
fprintf('Dataset generation completed!\n');
fprintf('Total samples: %d\n', nDatasetSize);
fprintf('Total time: %.2f min\n', tTotal/60);
fprintf('Average time per sample: %.2f sec\n', tTotal/nDatasetSize);
fprintf('Dataset saved to: %s\n', sDatasetDir);
fprintf('========================================\n');

% 데이터셋 요약 정보 저장
stDatasetInfo.nTotalSamples = nDatasetSize;
stDatasetInfo.nMinScatterers = nMinScatterers;
stDatasetInfo.nMaxScatterers = nMaxScatterers;
stDatasetInfo.nDepthRange_mm = [nMinDepth_mm, nMaxDepth_mm];
stDatasetInfo.nLateralRange_mm = [nMinLateral_mm, nMaxLateral_mm];
stDatasetInfo.ang_deg = ang_deg;
stDatasetInfo.nTxNum = nTxNum;
stDatasetInfo.generationTime = datetime('now');
stDatasetInfo.totalTimeSec = tTotal;

save([sDatasetDir 'dataset_info.mat'], 'stDatasetInfo');