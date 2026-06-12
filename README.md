# Reproducibility code for epithelial tissue necking

This repository contains the vertex-model simulation code, initial tissue
configuration, and visualization notebook used for the epithelial-tissue
necking simulations in the manuscript submitted to PLOS Computational Biology.

The code is a research fork of the vertex-model tutorial codebase by Rastko
Sknepnek:

- Upstream project: <https://github.com/sknepneklab/VMTutorial>
- License: MIT License, retained in `LICENSE`

The original `VMToolkit` package structure is retained. The paper-specific
workflow is organized in `TissueStretch/`.

## Repository Contents

```text
.
|-- VMToolkit/
|   |-- VM/                         C++ vertex-model engine and pybind11 module
|   |   `-- src/                    model forces, integrators, topology, dump code
|   `-- VMAnalysis/                 Python analysis helpers inherited from VMToolkit
|-- TissueStretch/
|   |-- honeycomb.json              initial tissue configuration
|   |-- run.py                      main tissue-stretch simulation script
|   `-- view_vtp.ipynb              notebook for opening and visualizing .vtp output
|-- config_builder/
|   |-- 10x3/honeycomb.json         small initial tissue geometry
|   `-- 30x10/honeycomb.json        large initial tissue geometry
|-- CMakeLists.txt
|-- pyproject.toml
|-- setup.py
`-- LICENSE
```

Build products, Python caches, and generated simulation outputs are not needed
for reproduction and can be regenerated.

## Requirements

The code was developed for Linux with Python 3.8 or newer. A Linux environment
is recommended because the Python extension is built from C++ with
CMake/scikit-build.

Required system tools:

- C++ compiler with C++14 support
- CMake
- Ninja

Required Python packages:

- `setuptools`
- `scikit-build`
- `pybind11`
- `cmake`
- `ninja`
- `boost`
- `numpy`
- `pyvista`

The build may also compile or locate VTK through the inherited VMTutorial build
configuration. If a system VTK is already installed, setting `VTK_DIR` before
installation can speed up the build.

## Installation

Clone the repository and install the package from the repository root:

```bash
git clone https://github.com/Yuanhe0517/Necking-of-epithelial-tissues.git
cd Necking-of-epithelial-tissues
python -m pip install --upgrade pip
python -m pip install .
```

Check that the compiled module imports correctly:

```bash
python -c "from VMToolkit.VM import Tissue, System, Force, Integrate; print('VMToolkit import OK')"
```

After modifying C++ source files in `VMToolkit/VM/src/`, rebuild and reinstall
the extension:

```bash
python -m pip install --force-reinstall --no-cache-dir .
```

## Running the Tissue-Stretch Simulation

The reproduction workflow is in `TissueStretch/`.

Run a short test:

```bash
cd TissueStretch
python run.py --input honeycomb.json --dt 0.1 --nrun 20 --loading-increment 0.03
```

Main arguments:

- `--input`: initial configuration JSON, default `honeycomb.json`
- `--dt`: initial integration time step, default `0.1`
- `--nrun`: number of outer loading/relaxation blocks, default `4000`
- `--seed`: random seed; if omitted, the code uses the internal default
- `--loading-increment`: boundary displacement applied after each saved relaxed state, default `0.03`

The `--input` argument can also point to one of the supplied geometries in
`config_builder/`:

```bash
python run.py --input ../config_builder/10x3/honeycomb.json --dt 0.1 --nrun 20 --loading-increment 0.03
python run.py --input ../config_builder/30x10/honeycomb.json --dt 0.1 --nrun 20 --loading-increment 0.03
```

The mechanical parameters are defined near the top of `TissueStretch/run.py`:

- area stiffness `kappa = 1.0`
- perimeter stiffness `gamma = 0.16`
- perimeter parameter `lambda = 0.56`
- boundary tension parameter `boundary_tension = 0`

T1 transitions are controlled by:

```python
topology.set_params({'min_edge_len': 0.1, 'new_edge_len': 0.12})
```

## Numerical Procedure

Each outer iteration calls:

```python
simulation.run(n)
```

where `n = 4000` for the first six outer iterations and `n = 3000` thereafter.
This is the maximum number of relaxation steps attempted in one call.

Mechanical equilibrium is detected in `IntegratorBrownian::step()` using the
successive total-energy difference:

```text
|E_k - E_{k-1}| < 1e-7
```

When this criterion is met, `simulation.run(n)` returns `true`, and `run.py`
writes the relaxed configuration. The loading increment is not applied until the
next call to the integrator. Thus the saved `.vtp` files correspond to relaxed
states before the next boundary displacement is applied.

If equilibrium is not reached within the maximum number of relaxation steps,
the time step is reduced in `Simulation::run()`:

```text
dt <- 0.9 dt
```

and the next outer iteration continues with the reduced time step.

## Outputs

The simulation writes output files in `TissueStretch/` when launched from that
directory.

Typical outputs:

- `cells_step_00000000.vtp`, `cells_step_00000001.vtp`, ...: relaxed cell configurations
- `600 steps.json`, `1200 steps.json`, ...: selected restart/state files
- `final steps.json`: final simulation state
- `data.csv`: T1 counts appended when equilibrium is reached

The `.vtp` files are VTK XML PolyData files. They can be opened in ParaView or
loaded in Python with PyVista.

## Visualizing VTP Files

Use the notebook:

```text
TissueStretch/view_vtp.ipynb
```

or load a `.vtp` file directly with PyVista:

```python
import pyvista as pv

mesh = pv.read("cells_step_00000000.vtp")
mesh.plot(show_edges=True)
```

## Initial Configurations

The default initial configuration copied into the simulation directory is:

```text
TissueStretch/honeycomb.json
```

Two additional initial tissue geometries are provided:

```text
config_builder/10x3/honeycomb.json
config_builder/30x10/honeycomb.json
```

The `10x3` configuration is useful for quick checks. The `30x10` configuration
is the larger tissue geometry used for production-scale stretch simulations.
Both files can be passed directly to `TissueStretch/run.py` with `--input`.

## Reproducibility Notes

For exact reproduction, record:

- Git commit hash of this repository
- input configuration file
- `dt`
- `nrun`
- `loading_increment`
- T1 parameters `min_edge_len` and `new_edge_len`
- convergence criterion
- Python version
- CMake version
- compiler version

## Citation and Attribution

If using this repository, cite the accompanying manuscript and acknowledge that
the simulation engine is derived from VMTutorial by Rastko Sknepnek:

<https://github.com/sknepneklab/VMTutorial>

The upstream VMTutorial code is distributed under the MIT License. Paper-specific
modifications and the tissue-stretch workflow were added by Yuan He for the
epithelial tissue necking study.
