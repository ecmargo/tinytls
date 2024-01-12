use ark_ec::CurveGroup;
use ark_std::UniformRand;
use rand::{CryptoRng, RngCore};
use transcript::IOPTranscript;

use crate::pedersen::CommitmentKey;
use crate::prover::{aes_prove, AesCipherWitness, AesKeySchWitness};
use crate::sigma::SigmaProof;
use crate::traits::{LinProof, ProofResult, TinybearProof};
use crate::verifier::{aes_verify, AesCipherInstance, AeskeySchInstance};
use crate::{aes, registry, u8msm};

#[inline]
pub fn aes128_prove<G: CurveGroup>(
    transcript: &mut IOPTranscript<G::ScalarField>,
    ck: &CommitmentKey<G>,
    message: [u8; 16],
    message_blinder: G::ScalarField,
    key: &[u8; 16],
    key_blinder: G::ScalarField,
) -> TinybearProof<G, impl LinProof<G>> {
    let witness =
        AesCipherWitness::<G::ScalarField, 11, 4>::new(message, key, message_blinder, key_blinder);
    aes_prove::<G, SigmaProof<G>, 11>(transcript, ck, &witness)
}

#[inline]
pub fn aes128ks_prove<G: CurveGroup>(
    transcript: &mut IOPTranscript<G::ScalarField>,
    ck: &CommitmentKey<G>,
    key: [u8; 16],
    key_blinder: G::ScalarField,
) -> TinybearProof<G, impl LinProof<G>> {
    let witness = AesKeySchWitness::<G::ScalarField, 11, 4>::new(&key, &key_blinder);
    aes_prove::<G, SigmaProof<G>, 11>(transcript, ck, &witness)
}

#[inline]
pub fn aes128_verify<G: CurveGroup>(
    transcript: &mut IOPTranscript<G::ScalarField>,
    ck: &CommitmentKey<G>,
    message_commitment: &G,
    round_keys_commitment: &G,
    ctx: [u8; 16],
    proof: &TinybearProof<G, impl LinProof<G>>,
) -> ProofResult {
    let instance =
        AesCipherInstance::<G, 11, 4>::new(message_commitment, round_keys_commitment, ctx);
    aes_verify::<G, 11>(transcript, ck, &instance, proof)
}

#[inline]
pub fn aes128ks_verify<G: CurveGroup>(
    transcript: &mut IOPTranscript<G::ScalarField>,
    ck: &CommitmentKey<G>,
    round_keys_com: G,
    proof: &TinybearProof<G, impl LinProof<G>>,
) -> ProofResult {
    let instance = AeskeySchInstance::<G, 11, 4>::new(&round_keys_com);
    aes_verify::<G, 11>(transcript, ck, &instance, proof)
}

#[inline]
pub fn aes256_prove<G: CurveGroup>(
    transcript: &mut IOPTranscript<G::ScalarField>,
    ck: &CommitmentKey<G>,
    message: [u8; 16],
    message_blinder: G::ScalarField,
    key: &[u8; 32],
    key_blinder: G::ScalarField,
) -> TinybearProof<G, impl LinProof<G>> {
    let witness =
        AesCipherWitness::<G::ScalarField, 15, 8>::new(message, key, message_blinder, key_blinder);
    aes_prove::<G, SigmaProof<G>, 15>(transcript, ck, &witness)
}

#[inline]
pub fn aes256ks_prove<G: CurveGroup>(
    transcript: &mut IOPTranscript<G::ScalarField>,
    ck: &CommitmentKey<G>,
    key: [u8; 32],
    key_blinder: G::ScalarField,
) -> TinybearProof<G, impl LinProof<G>> {
    let witness = AesKeySchWitness::<G::ScalarField, 15, 8>::new(&key, &key_blinder);
    aes_prove::<G, SigmaProof<G>, 15>(transcript, ck, &witness)
}

pub fn aes256_verify<G>(
    transcript: &mut IOPTranscript<G::ScalarField>,
    ck: &CommitmentKey<G>,
    m_com: &G,
    rk_com: &G,
    ctx: [u8; 16],
    proof: &TinybearProof<G, impl LinProof<G>>,
) -> ProofResult
where
    G: CurveGroup,
{
    let instance = AesCipherInstance::<G, 15, 8>::new(m_com, rk_com, ctx);
    aes_verify::<G, 15>(transcript, ck, &instance, proof)
}

pub fn commit_message<G: CurveGroup, const R: usize>(
    csrng: &mut (impl CryptoRng + RngCore),
    ck: &CommitmentKey<G>,
    m: [u8; 16],
) -> (G, G::ScalarField) {
    let message_offset = registry::aes_offsets::<R>().message;
    let m = m.iter().flat_map(|x| [x & 0xf, x >> 4]).collect::<Vec<_>>();
    let message_blinder = G::ScalarField::rand(csrng);
    let message_commitment =
        crate::u8msm::u8msm::<G>(&ck.vec_G[message_offset * 2..], &m) + ck.H * message_blinder;

    (message_commitment, message_blinder)
}

pub fn commit_aes128_message<G: CurveGroup>(
    csrng: &mut (impl CryptoRng + RngCore),
    ck: &CommitmentKey<G>,
    m: [u8; 16],
) -> (G, G::ScalarField) {
    commit_message::<G, 11>(csrng, ck, m)
}

pub fn commit_aes256_message<G: CurveGroup>(
    csrng: &mut (impl CryptoRng + RngCore),
    ck: &CommitmentKey<G>,
    m: [u8; 16],
) -> (G, G::ScalarField) {
    commit_message::<G, 15>(csrng, ck, m)
}

pub fn commit_aes128_key<G: CurveGroup>(
    csrng: &mut (impl CryptoRng + RngCore),
    ck: &CommitmentKey<G>,
    key: &[u8; 16],
) -> (G, G::ScalarField) {
    commit_round_keys(csrng, ck, &aes::aes128_keyschedule(key))
}

pub fn commit_aes256_keys<G: CurveGroup>(
    csrng: &mut (impl CryptoRng + RngCore),
    ck: &CommitmentKey<G>,
    key: &[u8; 32],
) -> (G, G::ScalarField) {
    commit_round_keys(csrng, ck, &aes::aes256_keyschedule(key))
}

fn commit_round_keys<G: CurveGroup, const R: usize>(
    csrng: &mut (impl CryptoRng + RngCore),
    ck: &CommitmentKey<G>,
    round_keys: &[[u8; 16]; R],
) -> (G, G::ScalarField) {
    let kk = round_keys
        .iter()
        .flatten()
        .flat_map(|x| [x & 0xf, x >> 4])
        .collect::<Vec<_>>();

    let key_blinder = G::ScalarField::rand(csrng);
    let round_keys_offset = registry::aes_offsets::<R>().round_keys * 2;
    let round_keys_commitment =
        u8msm::u8msm::<G>(&ck.vec_G[round_keys_offset..], &kk) + ck.H * key_blinder;

    (round_keys_commitment, key_blinder)
}
