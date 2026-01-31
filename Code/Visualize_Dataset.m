% Visualize_Dataset.m
%
% 생성된 데이터셋 시각화 및 검증
%
% 2025-01-31

clear; close all;

%% 설정
sDatasetDir = '..\dataset\';

%% 데이터셋 정보 로드
load([sDatasetDir 'dataset_info.mat']);
fprintf('Dataset Info:\n');
fprintf('  Total samples: %d\n', stDatasetInfo.nTotalSamples);
fprintf('  Scatterers: %d ~ %d\n', stDatasetInfo.nMinScatterers, stDatasetInfo.nMaxScatterers);
fprintf('  Depth range: %.1f ~ %.1f mm\n', stDatasetInfo.nDepthRange_mm(1), stDatasetInfo.nDepthRange_mm(2));
fprintf('  Lateral range: %.1f ~ %.1f mm\n', stDatasetInfo.nLateralRange_mm(1), stDatasetInfo.nLateralRange_mm(2));

%% 랜덤 샘플 시각화
nSamplesToShow = 5;
aIdx = randperm(stDatasetInfo.nTotalSamples, nSamplesToShow);

figure('Position', [100, 100, 1400, 800]);

for i = 1:nSamplesToShow
    dataIdx = aIdx(i);
    sDataID = sprintf('%06d', dataIdx);
    
    % 데이터 로드
    load([sDatasetDir 'bf_128ch\bf_' sDataID '.mat']);
    load([sDatasetDir 'bf_64ch\bf_' sDataID '.mat']);
    load([sDatasetDir 'metadata\meta_' sDataID '.mat']);
    
    % dB 변환
    mBf_128ch_dB = 20*log10(mBfData_128ch_norm + eps);
    mBf_64ch_dB = 20*log10(mBfData_64ch_norm + eps);
    
    % 128ch 이미지
    subplot(2, nSamplesToShow, i);
    imagesc(stMetadata.aX_128ch*1e3, stMetadata.aZ_128ch*1e3, mBf_128ch_dB);
    colormap(gray);
    clim([-60 0]);
    axis equal tight;
    title(sprintf('128ch #%d (%d pts)', dataIdx, stMetadata.nNumScatterers));
    xlabel('X (mm)');
    if i == 1
        ylabel('Z (mm)');
    end
    
    % 산란체 위치 표시
    hold on;
    plot(stMetadata.mScatPos_mm(:,1), stMetadata.mScatPos_mm(:,3), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    
    % 64ch 이미지
    subplot(2, nSamplesToShow, nSamplesToShow + i);
    imagesc(stMetadata.aX_64ch*1e3, stMetadata.aZ_64ch*1e3, mBf_64ch_dB);
    colormap(gray);
    clim([-60 0]);
    axis equal tight;
    title(sprintf('64ch #%d', dataIdx));
    xlabel('X (mm)');
    if i == 1
        ylabel('Z (mm)');
    end
    
    % 산란체 위치 표시
    hold on;
    plot(stMetadata.mScatPos_mm(:,1), stMetadata.mScatPos_mm(:,3), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    
end

sgtitle('Dataset Samples: 128ch (top) vs 64ch (bottom)', 'FontSize', 14);

%% 128ch vs 64ch 비교 (단일 샘플 상세)
figure('Position', [100, 100, 1200, 500]);

dataIdx = aIdx(1);
sDataID = sprintf('%06d', dataIdx);
load([sDatasetDir 'bf_128ch\bf_' sDataID '.mat']);
load([sDatasetDir 'bf_64ch\bf_' sDataID '.mat']);
load([sDatasetDir 'metadata\meta_' sDataID '.mat']);

mBf_128ch_dB = 20*log10(mBfData_128ch_norm + eps);
mBf_64ch_dB = 20*log10(mBfData_64ch_norm + eps);

subplot(1,3,1);
imagesc(stMetadata.aX_128ch*1e3, stMetadata.aZ_128ch*1e3, mBf_128ch_dB);
colormap(gray); clim([-60 0]); axis equal tight;
title('128ch (Ground Truth)');
xlabel('X (mm)'); ylabel('Z (mm)');
colorbar;

subplot(1,3,2);
imagesc(stMetadata.aX_64ch*1e3, stMetadata.aZ_64ch*1e3, mBf_64ch_dB);
colormap(gray); clim([-60 0]); axis equal tight;
title('64ch (Input)');
xlabel('X (mm)'); ylabel('Z (mm)');
colorbar;

% 차이 이미지 (리사이즈 필요할 수 있음)
subplot(1,3,3);
% 64ch를 128ch 크기로 리사이즈
mBf_64ch_resized = imresize(mBfData_64ch_norm, size(mBfData_128ch_norm));
mDiff = abs(mBfData_128ch_norm - mBf_64ch_resized);
imagesc(stMetadata.aX_128ch*1e3, stMetadata.aZ_128ch*1e3, mDiff);
colormap(hot); axis equal tight;
title('Difference |128ch - 64ch|');
xlabel('X (mm)'); ylabel('Z (mm)');
colorbar;

sgtitle(sprintf('Sample #%d: %d scatterers', dataIdx, stMetadata.nNumScatterers), 'FontSize', 14);

%% 데이터 크기 통계
fprintf('\n=== Data Size Statistics ===\n');
fprintf('128ch BF data size: [%d x %d]\n', size(mBfData_128ch_norm, 1), size(mBfData_128ch_norm, 2));
fprintf('64ch BF data size: [%d x %d]\n', size(mBfData_64ch_norm, 1), size(mBfData_64ch_norm, 2));

%% 저장
saveas(gcf, [sDatasetDir 'dataset_visualization.png']);
fprintf('\nVisualization saved to: %s\n', [sDatasetDir 'dataset_visualization.png']);
