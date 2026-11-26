//! The port against NTIA's C++ P.528 (`reference/ntia-p528`), compiled as this test's reference.
//!
//! A case agrees when both succeed with the same mode of propagation, the same warnings and every reported value
//! within tolerance, or both refuse it with the same error code.

mod compare;
mod edges;
mod examples;
mod grid;
