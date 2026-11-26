//! Compiles NTIA's C++ P.528 (`ntia/src`, `ntia/include`) into a static library for the reference bindings.
//!
//! The sources are compiled where they are, unmodified: their `#include "../../include/x.h"` paths resolve inside
//! `ntia/`. `DLLEXPORT` is `extern "C" __declspec(dllexport)`; on compilers that do not know `__declspec` it is defined
//! away, which leaves `extern "C"`.

use std::{fs, path::Path};

fn main() {
    println!("cargo:rerun-if-changed=ntia");

    let mut build = cc::Build::new();
    build.cpp(true).include("ntia/include").warnings(false);
    if !build.get_compiler().is_like_msvc() {
        build.define("__declspec(x)", Some(""));
    }
    for model in ["p528", "p676", "p835"] {
        let folder = Path::new("ntia/src").join(model);
        for entry in fs::read_dir(&folder).unwrap_or_else(|error| panic!("{}: {error}", folder.display())) {
            let path = entry.expect("an entry").path();
            if path.extension().is_some_and(|extension| extension == "cpp") {
                build.file(path);
            }
        }
    }
    build.compile("ntia_p528");
}
