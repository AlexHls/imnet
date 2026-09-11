# imnet

`imnet` is an interactive and headless tool for running nuclear reaction
networks, inspecting abundance evolution, and exporting reaction-flow state for
analysis. It uses Dear ImGui for the desktop interface and can be built against
one network backend at a time.

The network implementation is an external dependency. Yann and NuPPN sources
are not copied into this project; users provide paths to their local source
trees when configuring the build.

## Current Capabilities

- Build-time backend selection: `YANN` or `NUPPN`.
- GUI integration at fixed `rho`, `T`, and timestep.
- GUI integration to a final time, using backend substeps where available.
- GUI trajectory integration with cached rows that can be inspected like a
  trajectory.
- Nuclide-chart visualization, selected-isotope diagnostics, reaction flux
  arrows, trajectory plots, and abundance evolution plots.
- Headless single-step, final-time, and trajectory runs.
- CSV output for numerical runs.
- JSON state export containing species, trajectory/cache snapshots,
  abundances, energy generation, and reaction fluxes.

## Requirements

- CMake 3.16 or newer
- C compiler with C11 support
- C++ compiler with C++17 support
- GSL
- GLFW 3
- OpenGL
- pkg-config

NuPPN builds also require a Fortran compiler and the build tools expected by
the supplied NuPPN checkout.

On macOS with Homebrew:

```sh
brew install cmake pkg-config gsl glfw
```

On Debian/Ubuntu:

```sh
sudo apt-get install cmake pkg-config libgsl-dev libglfw3-dev gfortran
```

## Backend Source Trees

`imnet` expects backend source trees to be available locally:

- Yann: supply the directory containing `network.c`, `network.h`, and the
  associated Yann source/header files.
- NuPPN: supply the root NuPPN checkout. The build uses the `frames/ppn`
  source tree and the `frames/ppn/run_template` input directory by default.

Backend sources should be treated as external, read-only dependencies from the
application layer. Do not vendor proprietary backend code into public `imnet`
releases.

## Build

### Yann Backend

The Yann backend is the default backend. Supply the Yann source path explicitly:

```sh
cmake -S . -B build-yann \
  -DIMNET_NETWORK_BACKEND=YANN \
  -DIMNET_YANN_SOURCE_DIR=/path/to/yann
cmake --build build-yann
```

### NuPPN Backend

```sh
cmake -S . -B build-nuppn \
  -DIMNET_NETWORK_BACKEND=NUPPN \
  -DIMNET_NUPPN_SOURCE_DIR=/path/to/nuppn
cmake --build build-nuppn
```

The NuPPN build invokes NuPPN's own ppn makefile in the supplied source tree.
This may update generated files or objects inside that external checkout.

### Install

```sh
cmake --install build-yann --prefix /path/to/prefix
```

The executable is installed as:

```text
/path/to/prefix/bin/imnet
```

Example network data is installed under:

```text
/path/to/prefix/share/imnet/data
```

When running an installed binary, pass the desired data or NuPPN run directory
with `--data-dir`.

## Test

The current automated tests exercise the application layer with the Yann
backend:

```sh
cmake -S . -B build-yann \
  -DIMNET_YANN_SOURCE_DIR=/path/to/yann \
  -DBUILD_TESTING=ON
cmake --build build-yann
ctest --test-dir build-yann --output-on-failure
```

The NuPPN configuration currently has no registered CTest tests; use a build
smoke test plus representative headless runs.

## Run the GUI

From a source-tree build:

```sh
build-yann/src/imnet --data-dir data
```

For the NuPPN backend, `--data-dir` is the NuPPN ppn run directory. If omitted,
the compiled-in default is the supplied checkout's `frames/ppn/run_template`.

The GUI starts with these main windows:

- `Nuclide Chart`: isotope map, abundance coloring, and reaction flux arrows.
  Flux colors can be flat, linearly scaled, or logarithmically scaled.
- `Workflow`: conditions, composition editing, run mode, and cached-step
  browsing.
- `Isotope Information`: abundance, net rate, and reaction contributors for
  the selected isotope.

Common workflow:

1. Load or confirm network data from `File`.
2. Load an abundance file with `File -> Load Abundances`, or set a composition
   preset in `Workflow -> Run`.
3. Set `rho`, `T`, and `dt` in `Workflow -> Conditions`.
4. Choose `Fixed timestep`, `Final time`, or `Loaded trajectory` in
   `Workflow -> Run`.
5. Press `Run`.
6. Inspect cached steps in `Workflow -> Cached Steps`.

Useful views:

- `View -> Trajectory Plot` shows `rho(t)` and `T(t)`.
- `View -> Abundance Plot` plots selected isotope mass fractions over cached
  time. Isotopes can be added by search or from the selected nuclide. Log axes
  use a positive floor so zeros do not break plotting.
- `View -> Trajectory Editor` edits loaded `time`, `rho`, and `T` rows and can
  save a modified trajectory.
- `File -> Save State` writes a portable JSON analysis file.

For NuPPN builds, `View -> Settings` also exposes the ppn input text files used
by the run directory. Saving those inputs writes into the selected NuPPN run
directory.

## Headless Usage

Headless mode never opens an ImGui/GLFW window. It is intended for scripted
network comparisons and reproducible runs.

Show all options:

```sh
build-yann/src/imnet --help
```

### Single Timestep

```sh
build-yann/src/imnet --headless \
  --data-dir data \
  --abundances initial_abundances.txt \
  --rho 1e9 \
  --temp 5e9 \
  --dt 1e-6 \
  --output result.csv
```

### Final Time

```sh
build-yann/src/imnet --headless \
  --data-dir data \
  --abundances initial_abundances.txt \
  --rho 1e9 \
  --temp 5e9 \
  --dt 1e-6 \
  --final-time 1e-3 \
  --dt-max 1e-5 \
  --dt-factor 1.2 \
  --max-steps 10000 \
  --output results.csv
```

Final-time mode records backend history rows where available. These rows can be
used in the GUI or JSON state output to inspect substeps over time.

### Trajectory

```sh
build-yann/src/imnet --headless \
  --data-dir data \
  --abundances initial_abundances.txt \
  --trajectory data/trajectory_example.txt \
  --output results.csv
```

`--trajectory` and `--final-time` are separate run modes and cannot be used
together.

### Save JSON State

```sh
build-yann/src/imnet --headless \
  --data-dir data \
  --abundances initial_abundances.txt \
  --trajectory data/trajectory_example.txt \
  --save-state imnet_state.json
```

Headless mode accepts `--save-state` without `--output`. The saved JSON state
uses format tag `imnet-state-v1` and is intended for external analysis, not for
loading back into `imnet`.

## Input Files

### Abundances

Plain text, one species mass fraction per line:

```text
p 0.5
he4 0.5
```

Rules:

- Blank lines and lines beginning with `#` are ignored.
- Inline comments beginning with `#` are accepted.
- Each data row must contain `species mass_fraction`.
- Unknown species are reported as warnings.
- Missing species default to zero.
- Negative or non-finite abundance values are rejected.
- Thermodynamic keys such as `dt`, `rho`, and `temp` are rejected.

### Trajectories

Plain text columns:

```text
# time(s) rho(g/cm^3) temp(K)
0.0    1.0e9  5.0e9
1.0e-3 1.1e9  5.2e9
2.0e-3 1.2e9  5.4e9
```

Tracer-style files with metadata and `time T rho` columns are also accepted:

```text
# time T rho
AGEUNIT SEC
TUNIT T9K
RHOUNIT CGS
ID particle 46772
1.0e-3 5.0e-4 2.25e9
```

Rules:

- Blank lines and lines beginning with `#` are ignored.
- Each data row must contain either `time rho temp(K)` or `time T rho`.
- `rho` and `temp` must be positive.
- `TUNIT T9K`, `TUNIT T9`, and `TUNIT GK` values are converted to Kelvin.
- Time must be monotonic non-decreasing.
- Trajectory files do not contain species abundances.

## Output Files

Headless CSV output has this schema:

```text
step,time,rho,temp,dt,dedt,substeps,status,n,p,he4,...
```

The abundance columns follow the loaded backend species order. Numeric values
are written with high precision in scientific notation.

JSON state export contains:

- units
- selected step
- backend name and network file/run-directory paths
- species metadata and UI masks
- loaded trajectory rows
- cached integration steps
- abundance vectors
- energy generation rates
- directed reaction fluxes in `dY/dt`, `dX/dt`, and rate form

## CMake Options

| Option | Default | Meaning |
| --- | --- | --- |
| `IMNET_NETWORK_BACKEND` | `YANN` | Compile-time backend, either `YANN` or `NUPPN`. |
| `IMNET_YANN_SOURCE_DIR` | `./yann` | Path to a Yann source checkout. The default is for local development only. |
| `IMNET_NUPPN_SOURCE_DIR` | empty | Path to a NuPPN source checkout. Required for `NUPPN`. |
| `IMNET_NUPPN_RUN_DIR` | `frames/ppn/run_template` under the NuPPN checkout | NuPPN run directory used by default at runtime. |

## Caveats

- Backend choice is compile-time. Build separate binaries to compare Yann and
  NuPPN.
- Backend source trees are external dependencies. Public `imnet` releases
  should not include proprietary Yann or NuPPN source code.
- The application treats backend source as external, but the NuPPN build system
  itself may update generated files inside the supplied NuPPN checkout.
- NuPPN keeps process-global backend state; restart the application to switch
  to a different NuPPN run directory.
- Numerical results depend on backend version, rate files, species set,
  timestep controls, and input units. Record these with published runs.
- The GUI is an inspection tool, not a provenance system. Use headless commands
  and saved state files for reproducible analyses.
- JSON state files are analysis exports only; they are not loadable project
  files.

## Minimal Example

```sh
cat > /tmp/imnet_abundances.txt <<'EOF'
p 0.5
he4 0.5
EOF

build-yann/src/imnet --headless \
  --data-dir data \
  --abundances /tmp/imnet_abundances.txt \
  --trajectory data/trajectory_example.txt \
  --output /tmp/imnet_results.csv \
  --save-state /tmp/imnet_state.json
```
