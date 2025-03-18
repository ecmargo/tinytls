pub mod prover;
pub mod verifier;

// Re-export common functions
pub use crate::witness::trace::cipher::AesCipherTrace;
pub use crate::witness::trace::cipher::{aes128, aes256};
pub use crate::witness::trace::keyschedule::{aes128_keyschedule, aes256_keyschedule};
