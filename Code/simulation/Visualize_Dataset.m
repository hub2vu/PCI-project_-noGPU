% 1. 파일 로드 (Matlab File Exchange에서 'nrrdread' 함수 필요)
% annotation_25.nrrd 파일이 작업 폴더에 있어야 합니다.
data = nrrdread('annotation_25.nrrd');

% 2. 데이터 전처리 (uint32 -> double 변환)
% 메모리 절약을 위해 시뮬레이션할 2D 단면(Slice)만 추출합니다.
% 예: 뇌의 정중앙 단면(Sagittal View 기준)을 10개 층(250um 두께)만 추출
[dim1, dim2, dim3] = size(data);
slice_center = round(dim3 / 2); 
thickness = 10; 
sub_vol = double(data(:, :, slice_center : slice_center + thickness));

% 3. 산란체(Phantom) 생성
voxel_size = 25e-6; % 25 마이크로미터

% 조직이 있는 부분(0이 아닌 값)의 인덱스 추출
[x_idx, y_idx, z_idx] = ind2sub(size(sub_vol), find(sub_vol > 0));
region_ids = sub_vol(sub_vol > 0);

% 물리적 좌표계로 변환 (Meter 단위)
x_pos = (x_idx - mean(x_idx)) * voxel_size;
y_pos = (y_idx - mean(y_idx)) * voxel_size;
z_pos = (z_idx - mean(z_idx)) * voxel_size;

positions = [x_pos, y_pos, z_pos];

% 4. [핵심] 영역 ID -> 초음파 반사율(Amplitude) 매핑
% 각 고유 ID마다 랜덤한 반사 강도를 부여하여 대조도(Contrast) 생성
unique_ids = unique(region_ids);
num_ids = length(unique_ids);

% ID별 반사율 테이블 생성 (0~1 사이 랜덤 값)
reflectivity_map = containers.Map(unique_ids, rand(num_ids, 1));

% 각 산란체에 반사율 할당
amplitudes = zeros(size(region_ids));
for i = 1:length(region_ids)
    base_amp = values(reflectivity_map, {region_ids(i)});
    % 스페클(Speckle) 패턴을 만들기 위해 가우시안 노이즈 추가
    amplitudes(i) = base_amp{1} * (0.8 + 0.4 * randn()); 
end

% 5. Field II 시뮬레이션 실행 (예시)
% set_field('c', 1540);
% ... (Transducer 설정 코드) ...
% [rf_data, t0] = calc_scat(tx, rx, positions, amplitudes);