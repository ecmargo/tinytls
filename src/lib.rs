#![forbid(unsafe_code)]

/// Basic AES implementation.
pub mod aes;
/// Common AES functionalities.
pub mod aes_utils;
/// Basic AES implementation (key schedule).
pub mod aes_ks;
/// Basic AES implementation.
pub mod aes_plain;
/// Basic AES GCM implementation.
pub mod aes_gcm;
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
mod u8msm;
/// Verifier module.
mod verifier;

pub use exports::*;
