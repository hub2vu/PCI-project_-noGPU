function H = parse_nrrd_header_simple(headerText)
% parse_nrrd_header_simple
%  - NRRD header에서 필요한 필드만 파싱

    lines = splitlines(string(headerText));
    H = struct();
    H.dimension = NaN;
    H.sizes = [];
    H.type = "";
    H.encoding = "raw";

    for i = 1:numel(lines)
        s = strtrim(lines(i));
        if s == "" || startsWith(s, "#") || startsWith(s, "NRRD")
            continue;
        end

        parts = split(s, ":");
        if numel(parts) < 2, continue; end
        key = lower(strtrim(parts(1)));
        val = strtrim(join(parts(2:end), ":"));

        switch key
            case "dimension"
                H.dimension = double(str2double(val));
            case "sizes"
                nums = sscanf(val, '%d');
                H.sizes = nums(:)'; % row vector
            case "type"
                H.type = char(val);
            case "encoding"
                H.encoding = char(val);
        end
    end

    if isnan(H.dimension) || isempty(H.sizes) || strlength(H.type)==0
        error("NRRD header missing required fields (dimension/sizes/type).");
    end
end
