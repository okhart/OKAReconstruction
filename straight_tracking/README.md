# Track Reconstruction for Axion-like Particle Decays

## Overview

This project provides a complete track reconstruction framework for events from a detector simulation of axion-like particle decays. The input data consist of hits on several detector planes (Z positions with X and/or Y measurements). The goal is to reconstruct the original particle trajectories (tracks) in 3D space.

Two complementary algorithms are implemented:

- **Combinatorial Kalman Filter (CKF)** – a sequential, model-based approach that uses a Kalman filter to propagate the track state and perform hit association.
- **Hough Transform** – a global pattern recognition method that maps hits to lines in parameter space.

The code is designed to work both with **simulated data** (where truth information is available for validation) and **experimental data** (no truth, only reconstruction).

---

## Features

- Modular design: separate headers for track, Hough, Kalman, filtering, and configuration.
- Single unified `main.cpp`: algorithm selection via compile-time constants.
- Configurable parameters in `settings.h` (slope ranges, binning, χ² cuts, debug flags, etc.).
- Two data modes:
  - `Simulation` – enables comparison with true tracks and computes efficiency.
  - `Experimental` – runs reconstruction only, suitable for real data.
- Track filtering to remove duplicates and keep the best-quality tracks.
- Debug output and optional histogram saving for the Hough transform (controlled by `DEBUG`).
- Clear English comments and documentation.

---

## Repository Structure

```text
.
├── settings.h      # All configurable parameters
├── track.h         # Hit and track classes, line fitting
├── hough.h         # Hough transform implementation
├── kalman.h        # Combinatorial Kalman Filter
├── filter.h        # Track overlap removal
├── main.cpp        # Main driver program
└── README.md       # This file
```

---

## Requirements

- ROOT 6.x or later
- C++11 (or later) compiler (e.g. `g++`)
- Input ROOT file containing the required TTree structure

---

## Input Data Format

The program expects a ROOT file with a TTree named `decay_points`.

The tree must contain:

- **First 3 branches** – each is a `float[4]` containing the true track parameters `(a, b, c, d)` for one simulated particle. These branches are used only in `Simulation` mode.
- **Remaining branches** – each is a `TVector` (or compatible ROOT vector) with at least four elements:
  - `[0]` – Z coordinate
  - `[1]` – X coordinate (or NaN if unavailable)
  - `[2]` – Y coordinate (or NaN if unavailable)
  - `[3]` – measurement uncertainty (`sigma`)

Each remaining branch corresponds to one detector plane. The program automatically determines the number of detector planes from the tree.

---

## Output

The program prints reconstruction statistics to the console.

### Simulation mode

- Total processed events
- Events skipped because of:
  - ambiguous ZXY hits
  - insufficient number of reconstructed tracks
  - χ² cut
- Number of **GOOD EVENTS**, where every expected track is correctly matched to the corresponding true track
- Reconstruction efficiency:

```
GOOD EVENTS / TOTAL EVENTS
```

### Experimental mode

- Total processed events
- Events skipped because of:
  - ambiguous ZXY hits
  - insufficient number of reconstructed tracks
  - χ² cut
- Number of **GOOD EVENTS**, which are accepted ones
- Acceptance:

```
GOOD EVENTS / TOTAL EVENTS
```

The current implementation does not write reconstructed tracks to a ROOT file, although this can be added easily.

---

## Configuration (`settings.h`)

All configurable parameters are defined as constants inside `settings.h`.

### General

- `INPUT_FILE`
- `TREE_NAME`
- `MAX_EVENTS`

### Data Mode

- `DATA_MODE`
  - `Simulation`
  - `Experimental`

### Tracking Algorithm

- `TRACKING_ALGORITHM`
  - `Kalman`
  - `Hough`

### Debug

- `DEBUG`

### Combinatorial Kalman Filter

- `CKF_MIN_HITS`
- `CKF_CHI2_CUT`
- `CKF_PROCESS_NOISE`
- `CKF_AMBIGUITY_CHI2`
- `CKF_GREEDY`
- `CKF_GENEROUS`

### Hough Transform

- `HOUGH_A_MIN`
- `HOUGH_A_MAX`
- `HOUGH_B_MIN`
- `HOUGH_B_MAX`
- `HOUGH_A_BINS`
- `HOUGH_B_BINS`
- `HOUGH_OVERLAP`
- `HOUGH_THRESHOLD`

### Event Selection

- `CHI2_PER_NDF_CUT`
- `EXPECTED_N_TRACKS`

### Output

- `PRINT_STATISTICS`

---

## Reconstruction Pipeline

For every event the program performs the following steps:

1. Read an event from the ROOT tree.
2. Construct hit objects from detector branches.
3. Check for duplicated ZXY hits.
4. Run the selected tracking algorithm independently in the ZX and ZY projections.
5. Match ZX and ZY candidates using a common ZXY hit.
6. Build complete track candidates.
7. Remove duplicate tracks.
8. Apply χ² quality cuts.
9. Compare reconstructed tracks with truth information (Simulation mode only).
10. Print final statistics.

**Note:** the current implementation assumes exactly **five hits per track**:

- 2 ZX hits
- 2 ZY hits
- 1 ZXY hit

This assumption matches the detector geometry used in the simulation.

---

## Algorithms

### Combinatorial Kalman Filter

The CKF reconstructs tracks sequentially:

1. Build an initial seed from the first two detector planes.
2. Propagate the state to the next plane.
3. Calculate χ² for every compatible hit.
4. Select the best hit.
5. Update the state using the Kalman gain.
6. Continue until the last detector plane.

Two operating modes are available:

- **Greedy** – assigned hits are removed immediately.
- **Generous** – all compatible branches are explored.

---

### Hough Transform

The Hough transform represents every hit as a line in parameter space

```
coord = a · z + b
```

The algorithm:

- scans accumulator bins,
- collects all hits intersecting each bin,
- accepts bins with at least `HOUGH_THRESHOLD` hits,
- converts accepted bins into track candidates.

An overlap between neighboring bins is used to reduce edge effects.

---

### Track Fitting

Each reconstructed candidate is fitted independently in the ZX and ZY projections using weighted linear least squares.

The weights are

```
1 / sigma²
```

The fit returns:

- track parameters,
- covariance,
- χ².

---

### Duplicate Filtering

The `filterTracks()` function:

1. Sorts tracks by quality.
2. Prefers tracks with:
   - more hits,
   - smaller χ².
3. Removes tracks sharing too many hits with already accepted tracks.

---

## Compilation

Compile using ROOT with GNUmakefile on OKA server.

---

## Usage

Run the executable:

```bash
./track_reco
```

The program reads the input file specified in `settings.h` and prints reconstruction statistics.

To modify reconstruction parameters, edit `settings.h` and recompile.

---

## Future Work

- Save reconstructed tracks into a ROOT file.
- Support an arbitrary number of hits per track.
- Implement more advanced track quality estimators.
- Parallelize event processing.
- Include a more realistic detector response model.

---

## Cuts setup

- CHI2_PER_NDF_CUT and CKF_CHI2_CUT

Used according to statistical significance value.

For CKF_CHI2_CUT: ndf = 1.
For CHI2_PER_NDF_CUT (on the current OKA setup): ndf = 2.

| α                      | χ² (df=1) | χ² (df=2) | χ² / 2 (ndf=2) |
|------------------------|-----------|-----------|----------------|
| 0.50 (50%)             | 0.455     | 1.386     | 0.693          |
| 0.40 (40%)             | 0.708     | 1.833     | 0.917          |
| 0.333 (1/3)            | 0.954     | 2.197     | 1.099          |
| 0.30 (30%)             | 1.074     | 2.408     | 1.204          |
| 0.25 (25%)             | 1.323     | 2.773     | 1.386          |
| 0.20 (20%)             | 1.642     | 3.219     | 1.609          |
| 0.15 (15%)             | 2.072     | 3.794     | 1.897          |
| 0.10 (10%)             | 2.706     | 4.605     | 2.303          |
| 0.05 (5%)              | 3.841     | 5.991     | 2.996          |
| 0.025 (2.5%)           | 5.024     | 7.378     | 3.689          |
| 0.01 (1%)              | 6.635     | 9.210     | 4.605          |
| 0.001 (0.1%)           | 10.828    | 13.816    | 6.908          |

- CKF_AMBIGUITY_CHI2

Distribution is a difference of two χ². Calculations required.

## Noise

CKF_PROCESS_NOISE should be used to take into account multiple scattering.

In the current version it is set to zero, so noise was not taken into account.

## Instruction for debug with simulated data

1) Actual track parameters must be written in the first branches.

2) Hits must be written in the following branches in the same order (as track parameters).

## Recomended values of uncertancies (?)

1) DT - 0.2 sm
2) HODO:
2.1) Large pads - 1.36 sm
2.2) Small pads - 0.25 sm

---

## ! WARNING !

Hough-based algorithm don't work correctly in this version.

-

© 2026, Leonid Lapshin, lapshinleonid2005@gmail.com. All rights reserved.