%% 1. 데이터 로드 및 설정
filename = 'annotation_25.nrrd';
% nrrdread 함수가 필요합니다 (없으면 다운로드 혹은 addpath)
volData = double(nrrdread(filename));

% 해상도 설정 (Allen Mouse Brain Atlas는 25 마이크로미터)
voxel_size = 25e-6; 

%% 2. 관심 영역(ROI) 추출 - 중요!
% 전체를 다 하면 너무 느리므로, 뇌의 정중앙 단면(Coronal)을 
% 약 1mm 두께(40 복셀)만큼만 잘라냅니다.
[dimX, dimY, dimZ] = size(volData);

% Coronal View (앞뒤 방향) 기준으로 중앙 추출
slice_idx = round(dimY / 2); 
thickness = 40; % 40 voxels * 25um = 1mm 두께
slab_vol = volData(:, slice_idx-thickness/2 : slice_idx+thickness/2, :);

%% 3. 산란체 위치(Positions) 생성
% 값이 0이 아닌(뇌 조직이 있는) 복셀의 인덱스를 찾습니다.
[ix, iy, iz] = ind2sub(size(slab_vol), find(slab_vol > 0));
tissue_ids = slab_vol(slab_vol > 0); % 해당 위치의 조직 ID들

% 인덱스를 미터(meter) 단위로 변환
% Field II는 보통 Z축이 깊이 방향, X축이 측면 방향입니다.
% 데이터의 중심을 (0,0, Z_start)로 맞춥니다.

x_pos = (ix - mean(ix)) * voxel_size; 
y_pos = (iy - mean(iy)) * voxel_size; % 두께 방향 (Elevational)
z_pos = (iz - min(iz)) * voxel_size + 0.005; % 5mm 깊이부터 시작하도록 설정

% Field II 입력 형식: [N x 3] 행렬
positions = [x_pos, y_pos, z_pos];

%% 4. 산란체 강도(Amplitudes) 생성 - 핵심 기술
% 각 조직 ID마다 서로 다른 '평균 반사율'을 부여하고,
% 그 위에 '랜덤 스페클 노이즈'를 입힙니다.

unique_ids = unique(tissue_ids);
num_tissues = length(unique_ids);

% (1) ID별 기본 반사율 맵핑 (0.5 ~ 1.5 사이 랜덤 부여)
% 뇌실(Ventricle) 등 특정 부위를 0으로 만들고 싶다면 여기서 ID를 찾아 수정해야 함
reflectivity_map = containers.Map(unique_ids, 0.5 + rand(num_tissues, 1));

amplitudes = zeros(size(tissue_ids));

% (2) 각 산란체에 강도 부여
for i = 1:length(tissue_ids)
    % 해당 조직의 기본 반사율 가져오기
    base_amp = reflectivity_map(tissue_ids(i));
    
    % [중요] 가우시안 노이즈(randn)를 곱해 스페클 패턴 생성
    % 조직 내부는 균일하지 않고 미세한 산란체들이 무작위로 분포하기 때문
    amplitudes(i) = base_amp * (randn()); 
end

%% 5. 결과 확인 및 저장
fprintf('생성된 산란체 개수: %d 개\n', length(amplitudes));
fprintf('팬텀 크기: 가로 %.1f mm, 세로 %.1f mm, 깊이 %.1f mm\n', ...
    (max(x_pos)-min(x_pos))*1000, (max(y_pos)-min(y_pos))*1000, (max(z_pos)-min(z_pos))*1000);

% 3D 산란체 분포 시각화 (랜덤하게 10000개만 뽑아서 확인)
figure;
sample_idx = randsample(length(amplitudes), 1110000);
scatter3(positions(sample_idx,1)*1000, positions(sample_idx,2)*1000, positions(sample_idx,3)*1000, ...
    2, amplitudes(sample_idx), 'filled');
title('Field II Phantom Scatterer Distribution');
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
axis equal; view(3); colormap(gray); colorbar;

% 이제 positions와 amplitudes를 Field II의 calc_scat 함수에 넣으면 됩니다.
% 예: [rf_data, t0] = calc_scat(tx, rx, positions, amplitudes);