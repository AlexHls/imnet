# imnet Specification

## Purpose

`imnet` is an interactive nuclear network application built with Dear ImGui.
It uses an external nuclear network source tree as its computational core. The
backend source tree is treated as read-only application infrastructure and must
not be modified for this project.

The first implementation priority is correctness and failure resistance. User
friendliness is useful, but secondary to making the calculation path explicit,
recoverable, and usable without a GUI.

## Goals

1. Integrate the nuclear network for user-specified physical conditions and
   initial composition.
2. Visualize network results on the nuclide chart, including isotope abundance
   and production/destruction rate information.
3. Integrate to a final time at fixed `rho` and `T`, using `dt` as the initial
   timestep guess and bounded output timesteps.
4. Integrate a supplied thermodynamic trajectory `(time, rho, T)`, cache each
   step, and allow inspection of cached results.
5. Provide a headless mode that runs the same network logic from files without
   starting the GUI.

## Non-Goals

- Do not rewrite, fork, or patch backend source trees.
- Do not add a new GUI framework.
- Do not optimize for large production workflows before the basic calculation
  path is reliable.
- Do not hide recoverable errors. Show or print enough information to identify
  the invalid input or failed step.

## Existing Boundaries

- External backend source trees contain the network implementations and must
  remain unchanged.
- `src/network_wrapper.*` owns the C++ boundary around the selected backend.
- `src/app_state.*` owns application state, integration settings, and cached
  trajectory data.
- `src/file_io.*` owns simple text file parsing and writing.
- `src/ui_*` owns ImGui rendering and user interaction.
- `src/main.cpp` owns command-line parsing and dispatch between GUI and
  headless behavior.

## Core Data Model

### Species

Each isotope has:

- name
- proton number `Z`
- mass number `A`
- neutron number `N`
- abundance `X`
- optional UI state such as selected, active, and fixed

Species metadata is read from the initialized network backend, not from a
parallel source of truth.

### Physical Conditions

For direct integration:

- density `rho` in `g/cm^3`
- temperature `T` in `K`
- integration time step `dt` in `s`
- optional final integration time in `s`
- optional maximum final-time timestep in `s`
- optional final-time timestep growth factor
- optional maximum number of final-time integration steps
- abundance vector `X`

Inputs must be validated before integration:

- `rho > 0`
- `T > 0`
- `dt >= 0`
- final-time mode requires `final_time > 0`, `dt > 0`, `dt_max > 0`,
  `dt_factor > 0`, and `max_steps > 0`
- abundance vector size matches the network species count
- abundances are finite and non-negative

### Trajectory

A trajectory is an ordered table of:

```text
time rho temperature
```

The parser should accept comments and blank lines. Times must be finite and
monotonic non-decreasing. Each integration interval uses:

```text
dt_i = time[i + 1] - time[i]
rho_start = rho[i]
T_start = temperature[i]
rho_end = rho[i + 1]
T_end = temperature[i + 1]
X_i = cached composition at step i
```

The final trajectory row has no forward interval. It should still be available
for inspection, with `dt = 0` and no integration attempted from that final
point.

## Feature 1: Single Network Integration

### Required Behavior

- The user can specify `rho`, `T`, `dt`, and initial abundances.
- The application normalizes abundances only when explicitly requested or when
  the integration path documents that normalization is part of execution.
- Integration updates the abundance vector through the selected network wrapper.
- The UI displays the latest abundance values on the nuclide chart.
- The UI displays an isotope detail view for the selected or hovered isotope.
- The isotope detail view includes at minimum:
  - isotope name
  - abundance `X`
  - estimated net `dX/dt`
  - current `rho`, `T`, and `dt`

### Rate Information

The application must expose rate information that helps answer which reactions
produce or destroy a selected isotope.

Minimum acceptable first version:

- show net `dX/dt` for the isotope
- classify the sign as production, destruction, or near-zero

Preferred later version:

- list the top contributing reactions that produce the isotope
- list the top contributing reactions that destroy the isotope
- show each contribution under the current `rho`, `T`, and abundance vector
- rates and Q values need to be calculated from current `rho`, `T` 

## Feature 2: Trajectory Integration

### Required Behavior

- The user can load a thermodynamic trajectory file.
- The application validates the full trajectory before running it.
- The user can run the trajectory from the current initial abundance vector.
- The application caches one result per trajectory row.
- The user can step through cached trajectory rows in the UI.
- Selecting a cached step updates displayed `rho`, `T`, `dt`, abundance chart,
  isotope detail values, and energy generation information.
- The user can plot cached abundance evolution for chosen isotopes, with
  searchable isotope selection, optional logarithmic axes, and an optional
  legend.

### Cache Contents

Each cached trajectory step should store:

- trajectory index
- time
- rho
- temperature
- dt to next row
- abundance vector after integrating through that step's interval
- energy generation value returned by the integration
- backend substep count for the integration that produced the row
- integration status
- error message, if the step failed

### Failure Behavior

The default failure behavior should be stop-on-first-failure:

- keep all successfully cached prior steps
- mark the failed step with its error
- do not continue using a corrupted or unknown composition
- make the failure visible in GUI and headless output
- stop on first failure

## Feature 3: Headless Mode

### Required Behavior

Headless mode runs without constructing an ImGui window. It should support:

- loading network data files
- loading an initial abundance file
- loading a trajectory file
- running the trajectory
- writing results to an output file
- returning a non-zero exit code on validation or integration failure

Minimum command shape:

```text
imnet --headless \
  --data-dir ./data \
  --abundances initial_abundances.txt \
  --trajectory trajectory.txt \
  --output results.txt
```

For single-step headless integration, the same entry point may accept:

```text
imnet --headless \
  --data-dir ./data \
  --abundances initial_abundances.txt \
  --rho 1e6 \
  --temp 1e9 \
  --dt 1e-3 \
  --output result.txt
```

Final-time mode should accept the same `rho`, `temp`, and abundance inputs,
plus:

```text
imnet --headless \
  --data-dir ./data \
  --abundances initial_abundances.txt \
  --rho 1e6 \
  --temp 1e9 \
  --dt 1e-3 \
  --final-time 1 \
  --dt-max 1e-2 \
  --dt-factor 1.2 \
  --max-steps 10000 \
  --output result.txt
```

Final-time output should contain the backend-produced history rows needed for
inspection, with backend-reported substeps for each row.

Trajectory, final-time, and single-step mode should be viable in headless mode.
`--trajectory` and `--final-time` select separate multi-step modes and should
not be accepted together.

### Headless Output

The headless trajectory output should be CSV with a stable schema.
Minimum fields:

- step
- time
- rho
- temperature
- dt
- dedt
- substeps
- status
- abundance columns for every isotope

## File Formats

### Abundance Input

Existing simple format:

```text
h1 0.7
he4 0.3
```

Required parser behavior:

- comments and blank lines are ignored
- unknown species are warnings, not fatal errors
- missing species default to zero
- negative or non-finite abundances are fatal
- optional normalization is explicit

### Trajectory Input

Existing simple format:

```text
# time rho temp
0.0 1e6 1e9
0.1 9e5 9e8
```

Tracer-style format:

```text
# time T rho
AGEUNIT SEC
TUNIT T9K
RHOUNIT CGS
0.0 5.0e-4 1e9
```

Required parser behavior:

- comments and blank lines are ignored
- `AGEUNIT`, `TUNIT`, `RHOUNIT`, and `ID` metadata lines are accepted
- `time T rho` rows are accepted, with `T9/GK` converted to Kelvin
- malformed rows are fatal by default
- `rho <= 0` or `T <= 0` is fatal
- decreasing time is fatal

## GUI Requirements

The GUI should stay simple:

- one nuclide chart
- one integration/trajectory control panel
- one isotope information panel
- simple file dialogs or typed paths
- clear status text for success or failure

The nuclide chart should support:

- color by abundance
- indication of near-zero abundance
- selected isotope inspection
- optional indication of disconnected or inactive species

The UI should avoid complex workflow state. When an action fails, the previous
valid state remains visible.

## Implementation Priorities

1. Protect backend boundaries: all new behavior goes through C++ wrapper,
   app state, file I/O, CLI, or UI code.
2. Make headless dispatch real: `--headless` must not create a window.
3. Add explicit abundance input and result output options for headless mode.
4. Make trajectory validation strict and predictable.
5. Store trajectory cache records with status and error information.
6. Improve isotope rate diagnostics from net `dX/dt` toward per-reaction
  contributors if supported by the selected backend.
7. Keep GUI changes thin over shared app-state logic so GUI and headless modes
   use the same integration path.

## Acceptance Criteria

### Single-Step GUI

- Given valid network files and a valid abundance vector, clicking Integrate
  updates abundances without crashing.
- The selected isotope panel shows abundance and net rate data.
- Invalid `rho`, `T`, `dt`, or abundance input is rejected with a visible error.

### Trajectory GUI

- Given a valid trajectory, the application runs every interval and caches each
  row.
- The user can move between cached steps and see the corresponding abundance,
  `rho`, `T`, `dt`, and energy generation data.
- A failed trajectory step leaves prior cached steps inspectable.

### Headless

- `--headless` never attempts to open an ImGui/GLFW window.
- A valid trajectory plus abundance file produces an output file and exit code
  `0`.
- Invalid input produces no misleading success output and returns a non-zero
  exit code.
