# itu-p528: Recommendation ITU-R P.528 for Rust

A Rust implementation of **Recommendation ITU-R P.528-5**, *A propagation prediction method for aeronautical mobile
and radionavigation services using the VHF, UHF and SHF bands*. It gives the basic transmission loss between a terminal
on or near the ground and an airborne one, or between two airborne terminals, over a smooth earth, from 100 MHz to
30 GHz.

The crate ports the U.S. reference implementation by NTIA/ITS, the software ITU-R publishes with the Recommendation,
and returns the same answers: its tests compare every result with NTIA's C++ compiled from source.

## Usage

```rust
use p528::{P528Error, P528Input, Polarization, p528};

fn main() -> Result<(), P528Error> {
    let result = p528(&P528Input {
        distance_km: 100.0,
        h1_m: 100.0,
        h2_m: 15_000.0,
        frequency_mhz: 3_600.0,
        time_percent: 90.0,
        polarization: Polarization::Horizontal,
    })?;
    // 151.6 dB, one of NTIA's published examples.
    println!("{:.1} dB, {:?}", result.basic_transmission_loss_db, result.mode);
    Ok(())
}
```

For many distances between the same two heights, build a `P528Model` once and evaluate each distance with it. The
model does the work that depends only on the heights, frequency and polarisation once, and gives the same results as
`p528`:

```rust
use p528::{P528Error, P528Model, Polarization};

fn main() -> Result<(), P528Error> {
    let model = P528Model::create(20.0, 5_000.0, 3_000.0, Polarization::Vertical)?;
    for distance_km in [1.0, 10.0, 50.0] {
        match model.loss(distance_km, 50.0) {
            Ok(result) => println!("{distance_km} km: {:.1} dB", result.basic_transmission_loss_db),
            Err(error) => println!("{distance_km} km: {error}"),
        }
    }
    Ok(())
}
```

`P528Model` is `Send` and `Sync`. The model and `p528` keep no global state and contain no `unsafe` code.

### Features

| Feature | Default | Adds |
|---|---|---|
| `parallel` | No | `P528Model::par_losses`, which evaluates a slice of distances on rayon's thread pool. Adds the `rayon` dependency. |

Without features, the only dependency is `bitflags`.

## Inputs

The limits are the ones NTIA's reference enforces. Outside them, `p528`, `P528Model::create` and `P528Model::loss`
return a `P528Error` that carries NTIA's return code (`P528Error::ntia_code`). `p528` checks every input in NTIA's
order, so an input NTIA refuses gets NTIA's code, NaN and infinity included; a model checks the heights and frequency
when it is created, then the distance and time percentage.

| Field | Unit | Accepted | Notes |
|---|---|---|---|
| `distance_km` | km | 0 and above | Great-circle distance between the terminals. |
| `h1_m` | m | 1.5 to 80,000 | Low terminal. Must not be above `h2_m`. |
| `h2_m` | m | 1.5 to 80,000 | High terminal. |
| `frequency_mhz` | MHz | 100 to 30,000 | None. |
| `time_percent` | % | 1 to 99 | None. |
| `polarization` | None | `Horizontal`, `Vertical` | None. |

## Outputs

`P528Result` holds the basic transmission loss, the free-space loss and the gaseous absorption loss (all dB), the
distance the result is for, the elevation angle of the ray at the low terminal, the mode of propagation (line of sight,
diffraction or troposcatter) and `Warnings`, a set of NTIA's warning bits.

## Behaviour to know

The port computes what the reference computes, including where that is unexpected. Where the reference returns a value
that is not a loss, the port returns an error instead:

- **Long transhorizon paths are refused** with `P528Error::OutsideAtmosphere`. From about 2,400 km (2,400 km to
  3,500 km in the cases measured, depending on heights and frequency), the troposcatter common volume is above 100 km,
  the top of the P.835 reference atmosphere. The reference then uses its out-of-range codes as atmospheric values and
  returns NaN with a success code.
- **Terminals at the same height and 0 km apart are refused** with `P528Error::CoincidentTerminals`, NTIA's code 10.
  The reference defines that code but returns success with every loss 0.
- **NaN inputs, and an infinite distance, are refused** with `P528Error::NotFinite` when they pass NTIA's range checks.
  The reference accepts them and returns NaN.

Kept as the reference has it:

- **The terminals are not reordered.** A low terminal above the high one is refused (`P528Error::TerminalGeometry`,
  NTIA code 4). Pass the lower terminal as `h1_m`.
- **Heights above 20 km are accepted to 80 km** with `Warnings::HEIGHT_LIMIT_H1` or `HEIGHT_LIMIT_H2`. The in-force
  Recommendation covers 20 km; above it, the result is informative only.
- **A line-of-sight result is for the distance the ray geometry converged on**, which can differ from the input by up
  to 1 m (`P528Result::distance_km`).
- **One table interpolation reads a neighbouring entry less one.** In NTIA's `NakagamiRice`, the second interpolation
  between `K` rows uses `NakagamiRiceCurves[d_K - 1][d_p] - 1` where `[d_K - 1][d_p - 1]` would be the neighbouring
  entry. The port keeps it (`src/model/nakagami_rice.rs`).

## What the model includes

The gaseous absorption is computed as the reference computes it: a layered ray trace (Recommendation ITU-R P.676,
Annex 1) through the mean annual global reference atmosphere (Recommendation ITU-R P.835). Only the parts of P.676 and
P.835 that P.528 calls are ported. The model assumes a smooth earth; terrain is out of scope.

## Verification

`reference/ntia-p528` holds NTIA's C++ sources, unmodified, and compiles them with the `cc` crate into a test-only
dev-dependency. The test binary `ntia_reference` compares the port with it:

| Test | Cases |
|---|---|
| `examples::published_examples` | NTIA's four published examples, against the published values and against the reference. |
| `grid::generated_cases_agree` | Deterministic random cases, log-uniform in distance (0.01 km to 2,000 km), heights and frequency, uniform in time percentage, both polarisations. 2,000 by default; set `P528_CASES` for more. |
| `grid::tabulated_percentages_agree` | The 17 time percentages the tables are built on, 20 geometries each. |
| `edges::*` | Every refusal code, several failing checks at once, NaN and infinite inputs, coincident terminals, height warnings, distances around the line-of-sight limit, at zero and to 10,000 km, and the corners of the ranges. |
| `model::reused_model_agrees_with_reference` | One `P528Model` per height pair, evaluated at distances out of order and repeated, against a fresh reference call for each. |
| `model::parallel_batch_matches_sequential` | With `parallel`: `par_losses` returns the sequential results, in order. |

A case agrees when both succeed with the same mode and warnings and every loss within 0.01 dB, or both refuse it with
the same code, or the reference returns NaN or coincident-terminal zeros and the port the matching error. Measured on
Windows with the reference built by MSVC: 20,000 generated cases agree; over 3,000 cases the largest loss difference
is 3.4e-13 dB, and 97 % of the results are bit-identical.

Only the MSVC build of the reference has been compared. NTIA's code calls `abs` unqualified on `double` values, which
another C++ compiler and standard library could resolve to the integer overload; the comparison would then show it.

Not covered: no case found sets `Warnings::DIFFRACTION_TROPO_REGION` (a probe of 252 height, frequency and
polarisation combinations at 2,000 km found none), so the port's handling of that warning is ported but not compared.

Run the comparison (it needs a C++ compiler):

```sh
cargo test --release --all-features --test ntia_reference
```

## Performance

Measured on one Windows machine, release build, one thread.

| Call | Time |
|---|---|
| `p528`, 3,000 generated cases | 3.1 ms per call (NTIA's C++: 14 ms) |
| `P528Model::create`, low terminal 20 m, 3,000 MHz | 1.4 ms to 4.1 ms, rising with the high terminal's height |
| `P528Model::loss`, 20,000 distances from 0 km to 100 km, high terminal 3 km to 20 km (all line of sight) | 23 µs to 31 µs per distance |
| `P528Model::loss`, the same distances, high terminal 100 m (61 % line of sight) | 118 µs per distance |

The model computes once the terminals' geometry, the line-of-sight and transhorizon searches, and the gaseous
properties of every atmospheric layer between the two terminals. The line-of-sight absorption then traces each
distance's ray through the stored layers with the same arithmetic. A ray below the horizon at the low terminal, and
the transhorizon absorption, whose top height changes with distance, are still traced in full for each distance.

## Source

Ported from NTIA/ITS's repository, https://github.com/NTIA/p528, release v5.1 (commit
`98d5175766f3b6ce8d39dbb795a594aeab103393`), which implements P.528-5. The C++ and header files are identical at
commit `34c6e281ce9fe089ca06d915bb9061cb57750a52`, the latest when this port was made; only the licence wording,
citation and README changed after v5.1.

## Licence

The Rust code in this repository is under the MIT licence (`LICENSE`). The MIT licence does not cover
`reference/ntia-p528/ntia/`: those are NTIA's files, unmodified, under NTIA's terms in `reference/ntia-p528/LICENSE.md`.
NTIA's software is a work of the U.S. Government and is not subject to copyright protection in the United States.
`NOTICE` states how this port changes NTIA's software, as NTIA's terms ask.

A published package of this crate leaves out `reference/` and `tests/`: the tests need NTIA's C++ and a C++
compiler.
