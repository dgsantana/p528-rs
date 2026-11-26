# NTIA P.528 reference (test-only)

NTIA/ITS's C++ implementation of Recommendation ITU-R P.528-5, copied unchanged from https://github.com/NTIA/p528 at
release v5.1, commit `98d5175766f3b6ce8d39dbb795a594aeab103393` (`include/` and `src/` here as `ntia/include/` and
`ntia/src/`; `LICENSE.md` from commit `34c6e281ce9fe089ca06d915bb9061cb57750a52`, which corrects its wording). Works of
NTIA employees are not subject to copyright protection in the United States (Title 17 USC 105); NTIA's disclaimer is in
`LICENSE.md`.

`build.rs` compiles the sources where they are. `src/lib.rs` binds `P528_Ex` and mirrors its structures. This crate
is a dev-dependency only: the library's tests compare the Rust port with it. It is never built into the library.
