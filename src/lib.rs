#![forbid(unsafe_code)]

/// Basic AES GCM implementation.
pub mod aes_gcm;
/// Basic AES implementation (key schedule).
pub mod aes_ks;
/// Basic AES implementation.
pub mod aes_plain;
/// Common AES functionalities.
pub mod aes_utils;
/// AES keyschedule and cipher constraints.
mod constrain;
/// Interface functions publicly exposed.
mod exports;
/// Linear algebra toolkit.
mod linalg;
/// Core lookup sub-protocol.
mod lookup;
/// Pedersen commitment.
pub mod pedersen;
/// Prover module.
mod prover;
/// Helper module for the prover and verifier.
mod registry;
/// Core sigma protocols sub-protocols.
pub mod sigma;
/// Core sumcheck sub-protocol.
pub mod sumcheck;
/// Unit-tests.
#[cfg(test)]
mod tests;
/// Generic models used in the proof.
#[allow(non_snake_case)]
mod traits;
/// Fast MSM for u8 scalar elements.
mod umsm;
/// Verifier module.
mod verifier;
/// Witness gen module - GCM .
mod witness_gcm;
/// Witness gen module - KS .
mod witness_ks;
/// Witness gen module - AES .
mod witness_plain;

pub use exports::*;
