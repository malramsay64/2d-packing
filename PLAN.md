# Rewrite plan: 2D crystal packing optimiser in Rust (CLI + WASM)

## 1. What exists today

The search finds the densest periodic packing of a 2D shape by, for each
wallpaper group and each *isopointal* choice of Wyckoff sites, Monte Carlo
annealing the cell (a, b, γ) and the free site parameters (x, y, θ) while
rejecting overlapping configurations.

- `src/reference/*.c` is the original, working C implementation.
- `src/packing/*.cpp` is a partial C++/pybind11 port. It does not work yet.
- `src/pypacking/` holds Python wrappers and the Bilbao data scraper.

### Defects in the C++ port (fixed by the rewrite, not ported)

| Location | Problem |
|---|---|
| everywhere | `M_2_PI` is **2/π (0.637)**, not 2π. Angle ranges, modulo and incline maths are all wrong. |
| `initialise_structure` | `shared_ptr(&vector.back())` on a vector that reallocates: dangling pointers and a double free. `vector<Basis>` slices `CellLengthBasis`/`MirrorBasis`, and `get_random_value` is not virtual anyway. |
| `CellLengthBasis::get_random_value` | Returns `v + (1 + δ)` where it should return `v * (1 + δ)`. |
| MC loop | `load_basis(best)` runs inside the step loop because of a misplaced brace. `packing_max` is never updated, `count_replicas` stays 0, and `save_basis`/`load_basis` are undefined. Registering the spdlog `"console"` logger twice throws. |
| `check_intersection` | Skips pairs within the same Wyckoff orbit and a shape's own periodic images, so many overlaps are missed. |
| `intersects_with` | Position caches never add the shape centre. `prev` pointers are updated inside the inner loop, which pairs the wrong segments. |
| `check_for_intersection` | The shell test `2/π − γ < π/4` should compare against π − γ. A fixed 1–2 shells is a heuristic and not correct in general. |
| `generate_isopointal_groups` | `shape.mirrors % wyckoff.mirrors` divides by zero when the site has no mirrors. |
| random | A global unseeded `mt19937` mixed with `rand()`: not reproducible and not thread-safe. Shape flips (chirality) from the reference are missing. |

### Algorithmic limits inherited from the reference

1. **Radial shape representation.** A shape is r(θ) sampled at N angles, and
   orientations are snapped to that grid. Each overlap test crosses O(N²)
   segment pairs (N=360 → ~8k segment tests per pair) and is still approximate.
2. **Every pair, every image.** Each step re-checks every pair of shapes
   against a fixed set of periodic images. There is no symmetry reduction,
   no incremental update and no proper cutoff.
3. **Ad-hoc annealing.** The objective is a fixed geometric kT schedule on
   `1/φ`. Step sizes are fixed (`max_step_size=0.01`) and moves pick one
   random parameter. A fixed number of independent cycles run with no sharing.
4. **Uniform budget.** Every (group, isopointal) candidate gets the same
   step count, even though most of them are hopeless.
5. **Cell degeneracy.** Only the γ ∈ [π/4, 3π/4] clamp and an area-preserving
   rescale stop the cell collapsing. No lattice reduction is done.

## 2. Target architecture

A Cargo workspace:

```
crates/
  packing-core/   # pure library: no IO, no threads, no global state; builds for wasm32
  packing-cli/    # native binary: argument parsing, parallel search, file output
  packing-wasm/   # wasm-bindgen bindings + Web Worker entry point
  packing-data/   # wallpaper group / Wyckoff tables (generated, embedded)
web/              # minimal JS/TS front-end (worker pool, canvas/SVG viewer)
```

The core exposes a **resumable** optimiser
(`Optimizer::new(spec, seed)`, `opt.run(n_steps) -> Progress`,
`opt.best() -> Structure`). The CLI calls it in a loop. The browser calls it
in small slices so it can report progress and cancel without blocking.

### Libraries

| Purpose | Crate | Notes |
|---|---|---|
| Vector maths | `glam` (`DVec2`, `DMat2`) | Small and fast; compiles to wasm SIMD. |
| Deterministic trig | `libm` | Gives bit-identical `sin`/`cos` on native and wasm, so seeded runs reproduce across targets. |
| RNG | `rand_core` + `rand_xoshiro` (`Xoshiro256PlusPlus`) | Seedable, fast, no OS dependency. `getrandom` with the `js` feature is used only when no seed is given. |
| Serialisation | `serde`, `serde_json` | Structure I/O and the wasm boundary (with `serde-wasm-bindgen`). |
| Errors | `thiserror` | |
| CLI | `clap` (derive), `rayon`, `indicatif`, `tracing` + `tracing-subscriber` | Rayon runs candidates in parallel. |
| Output | `svg` crate | SVG renders of structures, matching the reference output. |
| WASM | `wasm-bindgen`, `js-sys`, `web-sys`, `console_error_panic_hook`; built with `wasm-pack` | Parallelism uses a Web Worker pool (one candidate per worker). An optional `wasm-bindgen-rayon` feature is available where COOP/COEP headers can be set. |
| Tests | `proptest`, `approx`, `insta` (snapshot), `wasm-bindgen-test` | |
| Benchmarks | `criterion` (native); a small timing harness under `wasm-bindgen-test` | |
| Overlap fallback (optional) | `parry2d-f64` | Only used as a test oracle or for non-convex shapes if the custom code is not enough. It is too large to ship in the hot path or the wasm bundle. |

The Wyckoff data stays sourced from the Bilbao server. The existing Python
scraper is reduced to a one-off generator that writes `wallpaper.json`, and
that file is committed. `packing-data` embeds it with `include_str!` and parses
it once with `serde`. A later option is a `build.rs` that turns it into static
tables.

## 3. Algorithmic changes

### 3.1 Shape representation and exact overlap

- Shapes are **vertex polygons**. A radial r(θ) input is still accepted and
  converted to vertices. Circles and other smooth shapes get a dedicated
  exact path: circle–circle and circle–polygon.
- **Convex shapes** are tested with the **Separating Axis Theorem**, which
  costs O(n+m) per pair (n=3–8 vertices for polygons), instead of O(N²) segment
  crossings.
- **Non-convex shapes** are split once into convex pieces (Hertel–Mehlhorn or
  ear-clipping plus merge). A pair overlaps if any piece pair overlaps, and
  bounding circles are checked per piece first.
- Orientation is **continuous**, so θ is no longer snapped to the sampling grid.
- Each test runs a cheap-to-exact cascade: circumcircle distance → inscribed
  circle overlap (an immediate *overlap*) → per-piece AABB → SAT.

### 3.2 Symmetry-reduced, incremental overlap checking

- Only the **asymmetric unit** is checked: one representative shape per
  occupied site against **all** images (every site × symmetry op × lattice
  translation, excluding only the identity image of itself). This is
  complete: it includes same-orbit images and self-images, which the C++
  port misses. It costs `k × M` tests instead of `M²` (k = number of sites,
  M = images).
- The set of lattice translations comes from the cutoff `|n·a + m·b| ≤ 2R_max`
  on the **reduced** cell (§3.3). This replaces the "1 or 2 shells" heuristic
  and is provably sufficient.
- Each image caches its transformed vertices. A move of site *i* recomputes
  only the images of orbit *i*, and only re-tests the pairs that involve
  orbit *i*. A cell move retests everything, with early exit on the first
  overlap. Most cell moves are compressions, and those are rejected cheaply.

### 3.3 Cell parametrisation and lattice reduction

- The cell is stored as an upper-triangular matrix `h = [[a, b cosγ], [0, b sinγ]]`,
  constrained by the group's lattice type: oblique, rectangular, square, or
  hexagonal (a=b, γ=60°).
- **Gauss (Lagrange) 2D lattice reduction** runs whenever the cell shear
  exceeds a threshold. It changes the basis, not the lattice, and site
  fractional coordinates are remapped with the same unimodular matrix.
  This removes the γ clamp and the area-preserving rescale hack.

### 3.4 Floppy-box isobaric Monte Carlo

This replaces the kT-on-1/φ schedule with a hard-particle **NPT
"floppy box" MC** (Filion et al. 2009; de Graaf et al. 2012), the
standard method for densest packings of hard polygons.

- Move types:
  - site translation (x, y in fractional coordinates);
  - site rotation;
  - discrete moves for mirror-constrained orientations (keeping the
    reference's 180° and 90° flips);
  - chirality flip;
  - cell moves: log-area scaling, shear, and aspect ratio.
- Acceptance for a cell move: `min(1, exp(−βP·ΔA + N·ln(A'/A)))`. Any
  overlap means rejection.
- **Pressure ramp** βP goes from low (fluid-like exploration) to very high
  (≈10⁶, effectively compression).
- **Adaptive step sizes**: each move type tunes its own step size during
  equilibration to hit about 30–40% acceptance, then freezes it so detailed
  balance holds in the production phase.
- **Final squeeze**: after the MC, a deterministic contact-finding stage
  bisects the uniform cell scale down to first contact, then runs
  coordinate-wise line searches on each free parameter at fixed overlap
  tolerance. This takes φ from "MC noise" to touching contacts.

### 3.5 Search over candidates

- **Correct enumeration**: isopointal combinations come from multisets of
  compatible Wyckoff sites, with a site-symmetry compatibility check that
  avoids the modulo-by-zero bug. They are deduplicated by an affine
  normaliser where the data supports it, otherwise by canonical letter
  sequence.
- **Analytical upper bound pruning**: φ ≤ N·A_shape / A_min_cell.
  A_min_cell comes from the requirement that images never come closer than
  `2R_in`. Candidates whose bound falls below the current best are skipped.
- **Successive halving** (racing): every candidate gets a short run. The top
  1/η of candidates by φ get η× more budget, and so on. This spends most of
  the compute on the few groups that matter.
- **Replicas per candidate**: several independent seeds per candidate. An
  optional **parallel tempering in pressure** mode swaps replicas between
  adjacent pressures. That replaces the reference's 32 isolated cycles.

## 4. Computational changes

- Hot loop: no heap allocation, no `Rc`/`Arc`, no virtual dispatch.
  State is plain arrays (struct-of-arrays for image vertices). Moves record
  an undo entry and revert it in place, instead of saving and restoring the
  whole basis.
- `f64` throughout. `libm` trig gives cross-target determinism, and wasm
  is built with `-C target-feature=+simd128` so glam can use SIMD.
- Release profile: `opt-level=3`, `lto="fat"`, `codegen-units=1`,
  `panic="abort"` for wasm, then `wasm-opt -O3`.
- Parallelism happens only at the candidate and replica level, which is
  embarrassingly parallel. The CLI uses `rayon`; the browser uses N Web
  Workers (`navigator.hardwareConcurrency`) pulling from a shared queue in JS.
  Each work item seeds its RNG as `hash(global_seed, group, isopointal, replica)`,
  so results do not depend on scheduling.
- The CLI writes JSON results (all candidates, ranked), SVG renders of the
  best N, and a CSV summary compatible with the reference columns.

## 5. Tests and verification

### 5.1 Unit and property tests (`cargo test`, `proptest`)

- **Geometry**:
  - SAT matches a brute-force oracle (segment crossing plus
    point-in-polygon, or `parry2d` in dev-deps) on random polygons;
  - overlap is symmetric, translation-invariant and rotation-invariant;
  - touching or ε-separated configurations are classified consistently;
  - convex decomposition reconstructs the original area.
- **Shapes**: area and inradius/circumradius of regular polygons match
  closed forms; radial→vertex conversion round-trips.
- **Wallpaper data**:
  - there are exactly 17 groups;
  - each group's operations are closed under composition modulo lattice
    translations;
  - Wyckoff multiplicities equal the orbit size of a generic point;
  - special positions are fixed by their site-symmetry group.
- **Cell**:
  - Gauss reduction preserves |det h| and the lattice (it is a unimodular
    integer transform);
  - remapped fractional coordinates give identical real positions;
  - image enumeration matches a brute-force 7×7 shell for random cells.
- **Overlap engine**:
  - the symmetry-reduced check equals a naive all-pairs check over an
    explicit 5×5 supercell;
  - the incremental check equals the full recheck after random move
    sequences;
  - move followed by revert restores state bit-for-bit.

### 5.2 Monte Carlo correctness

- The step-size controller converges to the target acceptance.
- Detailed-balance smoke test: hard disks at fixed low pressure reproduce
  the known equation of state (βP vs φ) within statistical error.
- Same seed gives an identical result: native vs native, and **native vs
  wasm** (run in CI through `wasm-bindgen-test` in headless Chromium,
  comparing a hashed trajectory).

### 5.3 Known-answer regression suite

The full search must reach each of these within tolerance (default 1e-3),
and every final structure must pass an **independent validator**: a naive
exact overlap check on a 5×5 supercell, with no shared code path with the
optimiser.

| Shape | Expected φ_max | Reference structure |
|---|---|---|
| Disk | π/(2√3) ≈ 0.9069 | hexagonal lattice |
| Square | 1.0 | p1 / p4 |
| Equilateral triangle | 1.0 | p2 (2 per cell) |
| Regular hexagon | 1.0 | p6 / p1 |
| Regular pentagon | (5−√5)/3 ≈ 0.9213 | p2 double lattice |
| Regular heptagon | ≈ 0.8926 | p2 double lattice |
| Ellipse/other | snapshot (`insta`) of best φ and group | guards against regressions |

### 5.4 Performance benchmarks (`criterion`, tracked in CI)

- Overlap pair test (ns/pair) for polygons of 3–8 vertices and for disks.
- MC steps/second for representative candidates (p2 pentagons, p6m disks).
- Wall-clock time to reach φ within 1e-3 of the known optimum for each
  regression shape, compared with the C reference built from `src/reference`.
- Acceptance targets:
  - ≥10× steps/s over the reference for polygons;
  - ≥5× lower time-to-solution;
  - wasm within 2× of native single-thread.

### 5.5 CI

`cargo fmt --check`, `cargo clippy -D warnings`, `cargo test --workspace`,
`wasm-pack test --headless --chrome crates/packing-wasm`, the regression
suite in release mode, and a benchmark smoke run. GitHub Actions replaces
Travis.

## 6. Delivery phases

0. **Baseline**: build `src/reference`, record φ and runtime for the
   regression shapes, and commit the numbers as `bench/baseline.json`.
1. **Data and core types**: export `wallpaper.json` and set up
   `packing-data`, the `Shape`, `Cell` and `Structure` types, and the data
   tests.
2. **Geometry**: SAT, convex decomposition, disk paths, and the geometry
   property tests.
3. **Overlap engine**: images, cutoff, reduction, incremental updates, and
   the equivalence tests.
4. **Optimiser**: floppy-box MC, adaptive steps, final squeeze, and the MC
   tests.
5. **Search and CLI**: enumeration, bounds, successive halving, rayon, and
   the outputs. The regression suite passes.
6. **WASM and web**: bindings, worker pool, viewer, and the
   cross-target determinism test.
7. **Clean-up**: remove the C++/pybind11 sources and submodules, and retire
   or slim the Python package. `src/reference` stays as a
   benchmark baseline until parity is shown.
