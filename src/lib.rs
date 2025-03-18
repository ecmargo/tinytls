#![forbid(unsafe_code)]

pub mod witness;

pub mod proof_system;

pub mod subprotocols;

/// Unit-tests.
#[cfg(test)]
mod tests;
/// Generic models used in the proof.
#[allow(non_snake_case)]
mod traits;
pub mod utils;

mod exports;
