# 20-Story 3D Frame OpenSeesPy Analysis

This folder contains an OpenSeesPy version of the 20-story 3D moment-resisting frame example in this repository.
The current default model uses a `3 x 3` bay plan.

For a Korean version of this document, see `readmekr.md`.

## Scripts

### `three_story_20_static.py`
Performs static and modal analysis.

- Loads:
  - Gravity: distributed floor load `-5 kN/m^2`
  - Lateral: top point load `-3000 kN` in `UX`
- Outputs:
  - `results/three_story_20_static_displacements.csv`
  - `results/three_story_20_static_reactions.csv`
  - `results/three_story_20_static_modal.csv`

### `three_story_20_dynamic.py`
Performs dynamic time-history analysis.

- Loads:
  - Gravity: distributed floor load `-5 kN/m^2`
  - No lateral point load
  - Ground motion: `0.3g` sine pulse
- Outputs:
  - `results/three_story_20_dynamic_top.csv`
  - `results/three_story_20_dynamic_story_drift.csv`

### `three_story_20_opensees.py`
Shared OpenSeesPy model builder used by the analysis scripts.

### `three_story_20_common.py`
Shared configuration constants and helpers:

- story and geometry defaults
- section property computation
- floor-load generation
- conversion from story count to modeled level count

### `visualize_structure.py`
Generates a matplotlib 3D rendering of the structure and saves `structure_view_3d.png`.

### `visualize_openseespy.py`
Uses `vfo` to show the OpenSeesPy model in an interactive viewer.

## Model Notes

- Geometry: `NUM_STORIES = 20`, which produces `21` modeled levels including the base.
- Bay layout: `NUM_BAYS_X = 3`, `NUM_BAYS_Z = 3`
- Story height: `4.0 m`
- Bay spacing: `4.0 m x 4.0 m`
- Total plan size: `12.0 m x 12.0 m`
- Section: `1.0 m x 1.0 m` square member
- Material: `E = 2.0e7 kPa`, `nu = 0.167`

This differs from the MATLAB `three_story_20` scripts, where `numFloors` is used as the number of modeled floor levels directly.
Bay counts can be changed in `three_story_20_common.py`.

## Usage

```bash
# Static analysis
python three_story_20_static.py

# Dynamic analysis
python three_story_20_dynamic.py
```

## Environment Setup

Create the environment from this folder:

```bash
cd /path/to/3D_Frame_FEM/openseespy_scripts
python -m venv .usr_venv
source .usr_venv/bin/activate   # Windows: .usr_venv\\Scripts\\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### Conda environment

```bash
cd /path/to/3D_Frame_FEM/openseespy_scripts
conda env create -f environment.yml
conda activate three-story-openseespy
```

Optional:

```bash
pip install vfo
```

## Windows (Conda): Tcl DLL Fix

OpenSeesPy on Windows may require `tcl86t.dll` and `tk86t.dll` under:

```text
<env>\\DLLs\\
```

If they only exist under:

```text
<env>\\Library\\bin\\
```

copy them after creating the environment:

```powershell
conda activate opspy39
conda install tk
copy $env:CONDA_PREFIX\\Library\\bin\\tcl86t.dll $env:CONDA_PREFIX\\DLLs\\
copy $env:CONDA_PREFIX\\Library\\bin\\tk86t.dll  $env:CONDA_PREFIX\\DLLs\\
python -c "import openseespy.opensees as ops; print('ok')"
```

The included `post_setup.ps1` automates that copy step.

## Recommended Python Version

- `3.9`: recommended
- `3.8`: generally stable
- `3.10+`: may need extra DLL fixes on Windows
