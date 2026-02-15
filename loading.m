function textSpinnerExample()
    % 텍스트 기반 회전 스피너 예시 (Command Window 전용)
    
    % 회전시키며 표시할 문자들
    progressSymbols = ['|','/','-','∖'];
    numSteps = 100;   % 반복 횟수(실제 작업 범위에 맞춰 조정)

    fprintf('작업 진행 중: ');
    for i = 1:numSteps
        % 실제 작업 (여기서는 예시로 0.1초 대기)


        % 현재 보여줄 스피너 문자 인덱스 계산
        symbolIdx = mod(i-1, length(progressSymbols)) + 1;
        
        % 스피너 문자 출력
        fprintf('%c', progressSymbols(symbolIdx));
        pause(0.1)
        % 한 글자 지우고 위치를 뒤로 이동(백스페이스)
        fprintf('\b');
    end
    fprintf('\n완료!\n');
end
