function [stScat, meta] = build_brain_phantom_fieldii(nrrd_annotation_path, out_mat_path, opts)
% build_brain_phantom_fieldii
% NRRD(annotation) -> brain mask -> Field II scatterer cloud 생성 + 캐시 저장
%
% 사용법:
%   [stScat, meta] = build_brain_phantom_fieldii(annPath, cachePath);        % 기본값
%   [stScat, meta] = build_brain_phantom_fieldii(annPath, cachePath, opts);  % 파라미터 지정
%
% opts (optional):
%   opts.Nscat     = 200000;
%   opts.vox_um    = 25;
%   opts.slab_mm   = 0.4;
%   opts.z0_mm     = 0.6;
%   opts.amp_scale = 1e26;
%   opts.amp_sigma = 0.6;
%   opts.seed      = 0;

    if nargin < 3 || isempty(opts)
        opts.Nscat     = 200000;
        opts.vox_um    = 25;
        opts.slab_mm   = 0.4;
        opts.z0_mm     = 0.6;
        opts.amp_scale = 1e26;
        opts.amp_sigma = 0.6;
        opts.seed      = 0;
    end

    % 재현성
    rng(opts.seed);

    % nrrdread_simple가 경로에 있어야 함
    V = double(nrrdread_simple(nrrd_annotation_path));
    mask = (V ~= 0);

    % slab 적용: y 중앙 근처만 남겨 2D 단면처럼 사용
    Ny = size(mask,2);
    y0 = round(Ny/2);

    vox = opts.vox_um * 1e-6; % meters
    slab_vox = max(1, round((opts.slab_mm*1e-3)/vox));

    ymin = max(1, y0 - slab_vox);
    ymax = min(Ny, y0 + slab_vox);
    mask(:, [1:ymin-1, ymax+1:end], :) = 0;

    idx = find(mask);
    if isempty(idx)
        error('Mask is empty after slab. slab_mm을 키우거나 nrrd가 올바른지 확인하세요.');
    end

    Npick = min(opts.Nscat, numel(idx));
    idx = idx(randperm(numel(idx), Npick));

    [ix, iy, iz] = ind2sub(size(mask), idx);

    Nx = size(mask,1);

    % voxel index -> meters (x/y 중앙정렬, z는 양수 depth로)
    x = (ix - (Nx+1)/2) * vox;
    y = (iy - (Ny+1)/2) * vox;
    z = (iz - min(iz)) * vox + opts.z0_mm*1e-3;

    stScat.mScatXYZPos = [x(:), y(:), z(:)];
    stScat.aScatMag    = opts.amp_scale * exp(opts.amp_sigma * randn(Npick,1));

    meta.mask     = mask;
    meta.opts     = opts;
    meta.vox_m    = vox;
    meta.size_vox = size(V);

    save(out_mat_path, 'stScat', 'meta', '-v7.3');
end
