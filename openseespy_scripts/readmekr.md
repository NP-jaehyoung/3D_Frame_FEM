# 20층 3D 프레임 OpenSeesPy 해석

이 폴더는 이 저장소의 20층 3D 모멘트 저항 프레임 예제를 OpenSeesPy로 구현한 버전입니다.
현재 기본 모델은 `3 x 3` 경간 평면을 사용합니다.

영문 문서는 `README.md`를 참고하세요.

## 스크립트

### `three_story_20_static.py`
정적 해석과 모드 해석을 수행합니다.

- 하중:
  - 중력 하중: 층별 분포하중 `-5 kN/m^2`
  - 수평 하중: 최상층 절점에 `UX` 방향 `-3000 kN` 점하중
- 결과 파일:
  - `results/three_story_20_static_displacements.csv`
  - `results/three_story_20_static_reactions.csv`
  - `results/three_story_20_static_modal.csv`

### `three_story_20_dynamic.py`
동적 시간이력 해석을 수행합니다.

- 하중:
  - 중력 하중: 층별 분포하중 `-5 kN/m^2`
  - 수평 점하중 없음
  - 지반운동: `0.3g` 사인 펄스
- 결과 파일:
  - `results/three_story_20_dynamic_top.csv`
  - `results/three_story_20_dynamic_story_drift.csv`

### `three_story_20_opensees.py`
해석 스크립트들이 공통으로 사용하는 OpenSeesPy 모델 생성 모듈입니다.

### `three_story_20_common.py`
공통 설정 상수와 보조 함수가 들어 있습니다.

- 층수 및 형상 기본값
- 단면 특성 계산
- 층하중 생성
- story 수를 실제 모델 level 수로 변환

### `visualize_structure.py`
구조물을 matplotlib 3D 그림으로 생성하고 `structure_view_3d.png`로 저장합니다.

### `visualize_openseespy.py`
`vfo`를 사용해 OpenSeesPy 모델을 인터랙티브 뷰어로 표시합니다.

## 모델 메모

- 형상: `NUM_STORIES = 20` 이며, base를 포함해 실제 모델 level은 `21`개입니다.
- 경간 구성: `NUM_BAYS_X = 3`, `NUM_BAYS_Z = 3`
- 층고: `4.0 m`
- 경간 간격: `4.0 m x 4.0 m`
- 전체 평면 크기: `12.0 m x 12.0 m`
- 단면: `1.0 m x 1.0 m` 정사각 단면
- 재료: `E = 2.0e7 kPa`, `nu = 0.167`

이 점은 MATLAB `three_story_20` 스크립트와 다릅니다. MATLAB 쪽은 `numFloors`를 실제 모델 floor level 수로 바로 사용합니다.
경간 수는 `three_story_20_common.py`에서 바꿀 수 있습니다.

## 사용 방법

```bash
# 정적 해석
python three_story_20_static.py

# 동적 해석
python three_story_20_dynamic.py
```

## 환경 설정

이 폴더에서 환경을 생성해서 사용하면 됩니다.

```bash
cd /path/to/3D_Frame_FEM/openseespy_scripts
python -m venv .usr_venv
source .usr_venv/bin/activate   # Windows: .usr_venv\\Scripts\\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### Conda 환경

```bash
cd /path/to/3D_Frame_FEM/openseespy_scripts
conda env create -f environment.yml
conda activate three-story-openseespy
```

선택 사항:

```bash
pip install vfo
```

## Windows (Conda): Tcl DLL 보정

Windows 환경에서는 OpenSeesPy가 아래 경로의 `tcl86t.dll`, `tk86t.dll`을 요구할 수 있습니다.

```text
<env>\\DLLs\\
```

만약 아래 경로에만 존재한다면:

```text
<env>\\Library\\bin\\
```

환경 생성 후 복사해주면 됩니다.

```powershell
conda activate opspy39
conda install tk
copy $env:CONDA_PREFIX\\Library\\bin\\tcl86t.dll $env:CONDA_PREFIX\\DLLs\\
copy $env:CONDA_PREFIX\\Library\\bin\\tk86t.dll  $env:CONDA_PREFIX\\DLLs\\
python -c "import openseespy.opensees as ops; print('ok')"
```

포함된 `post_setup.ps1` 스크립트로 이 복사 작업을 자동화할 수 있습니다.

## 권장 Python 버전

- `3.9`: 권장
- `3.8`: 비교적 안정적
- `3.10+`: Windows에서 추가 DLL 보정이 필요할 수 있음
