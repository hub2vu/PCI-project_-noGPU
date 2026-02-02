function stScat = make_brain_scat_from_annotation(nrrdPath, Nscat, slab_mm, z0_mm, amp_scale, amp_sigma)
% make_brain_scat_from_annotation
%  - annotation_25.nrrd에서 brain mask를 만들고 scatterer cloud 생성
%
% Inputs
%  nrrdPath  : 예) 'simulation\annotation_25.nrrd'
%  Nscat     : scatterer 개수 (권장: 5e4 ~ 3e5)
%  slab_mm   : y방향 slab 두께 (2D B-mode 느낌; 0.3~0.6mm 추천)
%  z0_mm     : 프로브면(z=0)과 뇌 윗면 사이 standoff (0.3~1.0mm 추천)
%  amp_scale : 산란 크기 스케일 (기본 1e26)
%  amp_sigma : 산란 랜덤 강도 (기본 0.6)
%
% Outputs
%  stScat.mScatXYZPos : [N x 3] (x,y,z) in meters
%  stScat.aScatMag    : [N x 1]

    if nargin < 5 || isempty(amp_scale), amp_scale = 1e26; end
    if nargin < 6 || isempty(amp_sigma), amp_sigma = 0.6;  end

    V = nrrdread_simple(nrrdPath);   % numeric 3D
    V = double(V);

    % annotation_25.nrrd는 25um spacing 가정
    vox_um = 25;
    vox = vox_um * 1e-6; % [m]

    % 1) brain mask (가장 안전한 방식: 0이 아닌 곳 = 뇌)
    mask = (V ~= 0);
    if ~any(mask(:))
        error("Mask empty: nrrd에 nonzero voxel이 없습니다. 파일/경로 확인");
    end

    % 2) y slab 제한 (2D처럼 보이게)
    Ny = size(V,2);
    y0 = round(Ny/2);
    slab_vox = max(1, round((slab_mm*1e-3)/vox));
    ymin = max(1, y0 - slab_vox);
    ymax = min(Ny, y0 + slab_vox);
    mask(:, [1:ymin-1, ymax+1:end], :) = 0;

    idx = find(mask);
    if isempty(idx)
        error("Mask empty after slab. slab_mm을 키우거나 파일 확인");
    end

    Npick = min(Nscat, numel(idx));
    idx = idx(randperm(numel(idx), Npick));
    [ix, iy, iz] = ind2sub(size(mask), idx);

    % 3) voxel index -> meters
    Nx = size(V,1);

    % 중심 정렬: x,y는 0 중심
    x = (ix - (Nx+1)/2) * vox;
    y = (iy - (Ny+1)/2) * vox;

    % z는 depth로: 최소를 0으로 맞춘 뒤 standoff 추가
    z = (iz - min(iz)) * vox + z0_mm*1e-3;

    stScat.mScatXYZPos = [x(:), y(:), z(:)];

    % 4) 산란 크기(텍스처): log-normal 느낌
    stScat.aScatMag = amp_scale * exp(amp_sigma * randn(Npick,1));
end
