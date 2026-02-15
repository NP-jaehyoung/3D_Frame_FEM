# 3D_Frame_FEM 저장소 요약 (2026-02-15 기준)

## 프로젝트 개요
이 저장소는 MATLAB 기반 2D/3D 구조물 FEM(유한요소법) 예제 및 문제 풀이 모음입니다.  
특히 3차원 프레임(보-기둥-거더) 모델링, 조립, 경계조건 적용, 하중해석, 고유진동 해석을 다루며, 최근에는 20층 강체 슬래브(디아프램) 조건이 반영된 예제가 추가되었습니다.

## 폴더/주요 파일

### 루트
- `ThreeDimFrame.m` : 기본 3D 프레임(샘플) 해석
- `ThreeDimFrame_2.m` : 소규모 3D 프레임 예제
- `ThreeDimFrame_12_18.m` : 12~18층 계열 변형 예제
- `ThreeDimFrame_20story.m` : 20층 3D 프레임 예제 (최신)
- `ThreeDimFrame_20story_ode45.m` : 지반 가속도 입력 기반 ODE45 동적응답 예제
- `drawingMesh.m`, `loading.m` : 메쉬/하중 관련 보조 함수
- `README.md` : 기존 짧은 기본 안내
- `.gitignore` : `.DS_Store` 등 불필요 파일 제외

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
- 결과를 CSV로 정리 저장(변위/반력/모드 진동수 테이블)
- 20층 모델(`ThreeDimFrame_20story.m`) 생성
- 강체 바닥(디아프램) 제약 적용
  - 각 층(2층~20층)에서 `UX`, `UZ`, `RY`를 동일화
  - 축소 좌표계 해석 후 전체 좌표계로 복원해 반력까지 계산
- 공통 모듈화 시작
  - `build3DFrameGeometry.m`, `assemble3DFrameMatrices.m`, `applyRigidDiaphragmConstraints.m`
  - `solveConstrainedStatic.m`, `solveEigenModes.m`, `ode45StateRhs3DFrame.m`
  - 두 메인 스크립트(`ThreeDimFrame_20story.m`, `ThreeDimFrame_20story_ode45.m`)가 동일한 함수를 공유

## 실행 가이드 (MATLAB)
1. MATLAB 작업 폴더를 저장소 루트로 이동
2. `ThreeDimFrame_20story.m` 실행
3. 결과는 `results/` 폴더의 CSV로 확인
   - `ThreeDimFrame_20story_displacements.csv`
   - `ThreeDimFrame_20story_reactions.csv`
   - `ThreeDimFrame_20story_modes.csv`
   - `ThreeDimFrame_20story_ode45_topHistory.csv`
   - `ThreeDimFrame_20story_ode45_storyDrift.csv`

## 라이선스
- `LICENSE` 참고

---
최종 확인일: 2026-02-15
