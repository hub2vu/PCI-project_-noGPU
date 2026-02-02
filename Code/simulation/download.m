function download_allen_data()
    % 시작 URL 설정
    baseUrl = 'https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/';
    
    % 현재 폴더를 다운로드 시작점으로 설정
    baseSavePath = pwd; 
    
    fprintf('=== Allen Brain Atlas 전체 다운로드 시작 (폴더 포함) ===\n');
    recursive_download(baseUrl, baseSavePath);
    fprintf('=== 모든 다운로드가 완료되었습니다 ===\n');
end

function recursive_download(currentUrl, currentSavePath)
    % 1. 해당 폴더가 없으면 생성
    if ~exist(currentSavePath, 'dir')
        mkdir(currentSavePath);
        fprintf('  [폴더 생성] %s\n', currentSavePath);
    end

    % 2. 웹페이지 내용 읽기
    try
        options = weboptions('Timeout', 60); % 타임아웃 넉넉하게
        htmlContent = webread(currentUrl, options);
    catch
        fprintf('  [접속 실패] %s (넘어갑니다)\n', currentUrl);
        return;
    end

    % 3. [파일] 찾아서 다운로드
    % 확장자가 nrrd, csv, json, txt 인 파일 찾기
    filePattern = 'href="([^"]+\.(nrrd|csv|txt|json))"';
    fileTokens = regexp(htmlContent, filePattern, 'tokens', 'ignorecase');
    
    if ~isempty(fileTokens)
        fileList = unique([fileTokens{:}]);
        for i = 1:length(fileList)
            fileName = fileList{i};
            cleanName = strrep(fileName, './', ''); % 경로 정리
            
            fileUrl = [currentUrl, cleanName];
            savePath = fullfile(currentSavePath, cleanName);
            
            % 파일이 이미 있으면 건너뛰기 (시간 절약)
            if exist(savePath, 'file')
                fprintf('    [Skip] 이미 있음: %s\n', cleanName);
            else
                fprintf('    [다운로드] %s ... ', cleanName);
                try
                    websave(savePath, fileUrl, options);
                    fprintf('완료\n');
                catch
                    fprintf('실패\n');
                end
            end
        end
    end

    % 4. [폴더] 찾아서 재귀 진입 (핵심 기능)
    % href가 / 로 끝나는 링크를 찾습니다.
    folderPattern = 'href="([^"]+/)"';
    folderTokens = regexp(htmlContent, folderPattern, 'tokens');
    
    if ~isempty(folderTokens)
        folderList = unique([folderTokens{:}]);
        
        for i = 1:length(folderList)
            subFolderName = folderList{i};
            
            % 무한 루프 방지: 부모 폴더(../)나 정렬 쿼리(?) 등은 무시
            if strcmp(subFolderName, '../') || startsWith(subFolderName, '/') || contains(subFolderName, '?')
                continue;
            end
            
            % 재귀 호출: 새 주소와 새 저장 위치를 만들어 다시 함수 실행
            newUrl = [currentUrl, subFolderName];
            newSavePath = fullfile(currentSavePath, subFolderName(1:end-1)); % 끝의 / 제거
            
            fprintf('  >> [폴더 진입] %s\n', subFolderName);
            recursive_download(newUrl, newSavePath); % <--- 자기 자신을 다시 호출
        end
    end
end