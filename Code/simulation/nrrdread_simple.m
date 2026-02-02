function V = nrrdread_simple(nrrdPath)
% nrrdread_simple (robust)
% - encoding: gzip/raw attached NRRD 지원
% - 헤더 끝(빈 줄)까지 fgetl로 읽고 ftell로 data 시작 위치를 정확히 잡음

    fid = fopen(nrrdPath, 'rb');
    if fid < 0
        error("Cannot open: %s", nrrdPath);
    end

    % ---- 1) Read header lines until blank line ----
    headerLines = {};
    while true
        tline = fgetl(fid);
        if ~ischar(tline)
            fclose(fid);
            error("NRRD header ended unexpectedly (EOF).");
        end
        if isempty(tline)
            break; % blank line => end of header
        end
        headerLines{end+1} = tline; %#ok<AGROW>
    end
    headerText = strjoin(headerLines, newline);

    H = parse_nrrd_header_simple(headerText);

    if H.dimension ~= 3
        fclose(fid);
        error("Only 3D supported. dimension=%d", H.dimension);
    end

    % ---- 2) Read remaining bytes (data) ----
    dataBytes = fread(fid, inf, '*uint8');
    fclose(fid);

    dtype = nrrd_type_to_matlab_class(H.type);
    nElem = prod(H.sizes);

    enc = lower(strtrim(H.encoding));

    if contains(enc, 'gzip')
        % dataBytes가 gzip 스트림이어야 함 (매직 1F 8B)
        if numel(dataBytes) < 2 || ~(dataBytes(1)==hex2dec('1F') && dataBytes(2)==hex2dec('8B'))
            error("Data after header is not gzip stream. (magic 1F 8B not found)");
        end

        % gzip bytes -> temp.gz -> gunzip -> read
        tmpGz = [tempname, '.gz'];
        fid = fopen(tmpGz, 'wb');
        fwrite(fid, dataBytes, 'uint8');
        fclose(fid);

        outDir = fileparts(tmpGz);
        gunzip(tmpGz, outDir);

        % gunzip 결과 파일 찾기
        d = dir(outDir);
        d = d(~[d.isdir]);
        [~, k] = max([d.datenum]);
        outPath = fullfile(d(k).folder, d(k).name);

        fid = fopen(outPath, 'rb');
        Vflat = fread(fid, nElem, ['*' dtype]);
        fclose(fid);

        % cleanup
        try delete(tmpGz); end %#ok<TRYNC>
        try delete(outPath); end %#ok<TRYNC>

    elseif contains(enc, 'raw')
        % raw attached
        % typecast는 uint8 row vector가 필요
        dataBytes = uint8(dataBytes(:));
        Vflat = typecast(dataBytes, dtype);
        Vflat = Vflat(:);

        if numel(Vflat) < nElem
            error("Not enough raw data: expected=%d got=%d", nElem, numel(Vflat));
        end
        Vflat = Vflat(1:nElem);

    else
        error("Unsupported encoding: %s", H.encoding);
    end

    V = reshape(Vflat, H.sizes);
end

function cls = nrrd_type_to_matlab_class(t)
    t = lower(strtrim(t));
    if any(strcmp(t, {'uchar','uint8','unsigned char'})), cls='uint8'; return; end
    if any(strcmp(t, {'char','int8','signed char'})), cls='int8'; return; end
    if any(strcmp(t, {'short','int16'})), cls='int16'; return; end
    if any(strcmp(t, {'ushort','uint16','unsigned short'})), cls='uint16'; return; end
    if any(strcmp(t, {'int','int32'})), cls='int32'; return; end
    if any(strcmp(t, {'uint','uint32','unsigned int'})), cls='uint32'; return; end
    if any(strcmp(t, {'float','single'})), cls='single'; return; end
    if any(strcmp(t, {'double'})), cls='double'; return; end
    error("Unsupported NRRD type: %s", t);
end
