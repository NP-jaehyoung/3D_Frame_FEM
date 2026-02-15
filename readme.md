# 3D_Frame_FEM 저장소 요약 (2026-02-15 기준)

## 프로젝트 개요
MATLAB 기반 2D/3D 구조물 FEM(유한요소법) 예제 및 문제 풀이 모음입니다.  
3차원 프레임(보-기둥-거더) 모델링, 조립, 경계조건 적용, 하중해석, 고유진동 해석을 다룹니다.

최근 `three_story_20/ThreeDimFrame_20story.m`와 `three_story_20/ThreeDimFrame_20story_ode45.m`를 중심으로 구조화되어 있으며, 20층 고층 모델에서 강체 바닥(디아프램), 층별 분포하중, 동적응답(ODE45) 예제를 함께 제공합니다.

## 폴더/주요 파일

### 루트
- `ThreeDimFrame.m` : 기본 3D 프레임(샘플) 해석
- `ThreeDimFrame_2.m` : 소규모 3D 프레임 예제
- `ThreeDimFrame_12_18.m` : 12~18층 계열 변형 예제
- `drawingMesh.m`, `loading.m` : 메쉬/하중 관련 보조 함수
- `readme.md` : 기존 짧은 기본 안내
- `.gitignore` : `.DS_Store` 등 불필요 파일 제외

### `three_story_20/`
- `ThreeDimFrame_20story.m` : 20층 3D 프레임 정적 해석(변위/반력/고유진동수)
- `ThreeDimFrame_20story_ode45.m` : 지반 가속도 입력 기반 ODE45 동적응답 예제
- `utils/` : 20층 예제에서 사용하는 공통 유틸리티 함수
  - `build3DFrameGeometry.m`
  - `assemble3DFrameMatrices.m`
    - `floorLoads`(층별 하중 벡터)와 `floorLoadDof`(예: `"UX"`, `"UY"`, `"RZ"`)를 지원합니다.
    - `topLoad`(기존 상부 집중하중) + 층별 하중 동시 조합이 가능합니다.
  - `applyRigidDiaphragmConstraints.m`
    - 층별 바닥의 강체 조건을 `UX, UY, UZ, RX, RY, RZ`로 묶어 회전(비틀림) 변형까지 반영합니다.
  - `solveConstrainedStatic.m`
  - `solveEigenModes.m`
  - `ode45StateRhs3DFrame.m`

### `Problem/`
교재/연습형 문제 스크립트 모음
- `problem1.m` ~ `problem17a.m` (vibration, buckling, FGM, plate/truss/frame 등 범위별 실험)
- 주요 축
  - 정적/동적 해석 문제군: `problem3vib.m`, `problem5vib.m`, `problem7vib.m`, `problem9vib.m`, `problem11vib.m`, `problem16vibrations*.m`
  - 좌굴/안정성: `problem16Buckling.m`, `problem9buk.m`
  - FGM/재료 모델: `problem16fgm*.m`, `problem16timeReddy.m`
  - 판/격자/기타 FEM 유틸리티 기반 문제들: `problem16*.m` 계열

### `matlab_codes_fem_book/`
재사용 가능한 FEM 함수 라이브러리
- 강성/질량 조립
  - `formStiffness3Dframe.m`, `formStiffness2Dframe.m`, `formMass3Dframe.m`, `formMass2Dframe.m`
- 빔/트러스/판/Mindlin 요소 관련
  - `formStiffnessTimoshenkoBeam.m`, `formStiffnessBernoulliBeam.m`, `formStiffness3Dtruss.m`, `formStiffness2Dtruss.m`
  - `formMassMatrixMindlin*.m`, `formStiffnessMatrixMindlin*.m`, `formMass3Dtruss.m` 등
- 바운더리/해석 보조
  - `EssentialBC.m`, `solution.m`, `eigenvalue.m`, `outputDisplacementsReactions.m`
- 형상/적분/재료
  - `Jacobian*.m`, `gaussQuadrature.m`, `shapeFunction*.m`, `rectangularMesh.m`
  - `srinivasMaterial.m`, `reddy...`, `liewMaterial.m` 등

## 최근 변경 포인트 (요약)
- `.DS_Store` 무시 규칙 적용 (`.gitignore`)
- 결과를 CSV로 정리 저장(변위/반력/모드 진동수/외력/회전/층하중 테이블)
- 20층 모델(`ThreeDimFrame_20story.m`) 생성
- 강체 바닥(디아프램) 제약 적용
  - 각 층(2층~최상층)에서 바닥 평면 강체 변형 관계 적용  
    - `ux_i = uxc - dz*Ry + dy*Rz`
    - `uy_i = uyc + dz*Rx - dx*Rz`
    - `uz_i = ucz - dy*Rx + dx*Ry`
  - 축소 좌표계 해석 후 전체 좌표계로 복원해 반력까지 계산
- 공통 모듈화 시작
  - `three_story_20/utils/` 폴더의 스크립트를 통해 `ThreeDimFrame_20story.m`와
    `ThreeDimFrame_20story_ode45.m`가 공통 함수를 공유
- 단면관성 자동 계산: `Iy`, `Iz`는 `width`, `depth`에서,  
  `J`는 사각단면 근사식으로 `width/depth`에 따라 자동 산정

## 20층 예제 사용 포인트
- 기본 형상/해석 옵션
  - 층수: `numFloors = 20`
  - 층고: `floorHeight = 4.0` (m)
  - 평면 치수: `spanX = 4.0`, `spanZ = 4.0` (m)
- 지반 가속도 입력은 `groundAcceleration` 벡터로 설정 후 `ode45`에 전달
- 층하중 사용 예시
  - 모든 층 동일 하중: `floorLoads = -3000 * ones(numFloors,1)`
  - 층별 가중: `floorLoads = [0; -1000; -1500; ...]`
  - 방향 지정: `assemble3DFrameMatrices(..., floorLoads, "UX")` (기본값 UX)
- 출력 CSV (실행 시 `three_story_20/results/`)
  - `ThreeDimFrame_20story_displacements.csv`
  - `ThreeDimFrame_20story_reactions.csv`
  - `ThreeDimFrame_20story_modes.csv`
  - `ThreeDimFrame_20story_floorLoads.csv`
  - `ThreeDimFrame_20story_external_loads.csv`
  - `ThreeDimFrame_20story_node_rotations.csv`
  - `ThreeDimFrame_20story_ode45_topHistory.csv`
  - `ThreeDimFrame_20story_ode45_storyDrift.csv`

## 실행 가이드 (MATLAB)
1. MATLAB 작업 폴더를 `three_story_20` 폴더로 이동
   - `cd(".../3D_Frame_FEM/three_story_20")`
2. 정적 해석: `ThreeDimFrame_20story`
3. 지반 가속도 동적해석: `ThreeDimFrame_20story_ode45`
4. 생성 파일은 `three_story_20/results/`에서 확인 (위 표 참조)

## 라이선스
- `LICENSE` 참고

---
최종 확인일: 2026-02-15
