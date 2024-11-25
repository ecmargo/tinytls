use std::process::Output;

use crate::linalg::{self, SparseMatrix};
use crate::witness;
use crate::witness::trace::utils::{rotate_right_inplace, RJ2};
use crate::witness::{registry, trace::utils};
use ark_ff::{Field, Zero};

pub fn aes_trace_to_needles<F: Field, const R: usize>(
    output: &[u8; 16],
    src: &[F],
    [c_xor, c_xor2, c_sbox, c_rj2]: [F; 4],
) -> (Vec<F>, F) {
    let reg = registry::aes_offsets::<R>();

    let mut dst = vec![F::zero(); reg.witness_len * 2];
    let mut offset = 0;
    cipher_sbox::<F, R>(&mut dst, src, c_sbox);
    offset += 16 * (R - 1);
    cipher_rj2::<F, R>(&mut dst, &src[offset..], c_rj2);
    offset += 16 * (R - 2);
    cipher_mcol::<F, R>(&mut dst, &src[offset..], c_xor, c_xor2);
    offset += 16 * (R - 2) * 4 * 2;
    let constant_term = cipher_addroundkey::<F, R>(output, &mut dst, &src[offset..], c_xor, c_xor2);

    (dst, constant_term)
}

pub fn aes_keysch_trace_to_needles<F: Field, const R: usize, const N: usize>(
    src: &[F],
    [c_xor, c_xor2, c_sbox, _c_rj2]: [F; 4],
) -> (Vec<F>, F) {
    let registry = registry::aes_keysch_offsets::<R, N>();
    let mut dst = vec![F::zero(); registry.witness_len * 2];
    let mut offset: usize = 0;
    crate::constrain::ks_lin_sbox_map::<F, R, N>(&mut dst, src, c_sbox);
    offset += 4 * (R - N / 4);
    let constant_term =
        crate::constrain::ks_lin_xor_map::<F, R, N>(&mut dst, &src[offset..], [c_xor, c_xor2]);
    (dst, constant_term)
}

pub fn _aes_gcm_block_trace_to_needles<F: Field, const R: usize>(
    aes_output: &[u8; 16],
    _final_xor: &[u8; 16],
    src: &[F],
    [c_xor, c_xor2, c_sbox, c_rj2]: [F; 4],
) -> (Vec<F>, F) {
    let regions = registry::aes_gcm_block_offsets::<R>();

    let mut dst = vec![F::zero(); regions.witness_len * 2];
    let mut offset = 0;
    cipher_sbox::<F, R>(&mut dst, src, c_sbox);
    offset += 16 * (R - 1);
    cipher_rj2::<F, R>(&mut dst, &src[offset..], c_rj2);
    offset += 16 * (R - 2);
    cipher_mcol::<F, R>(&mut dst, &src[offset..], c_xor, c_xor2);
    offset += 16 * (R - 2) * 4 * 2;
    let constant_term =
        cipher_addroundkey::<F, R>(aes_output, &mut dst, &src[offset..], c_xor, c_xor2);
    gcm_final_xor::<F, R>(&mut dst, &src[offset..], c_xor, c_xor2);
    //not including the constant term here because im not sure its needed?
    (dst, constant_term)
}

pub fn cipher_mcol<F: Field, const R: usize>(dst: &mut [F], v: &[F], r: F, r2: F) {
    let identity = (0..16).collect::<Vec<_>>();
    let registry = registry::aes_offsets::<R>();

    let mut aux_m_col = vec![identity; 4];
    utils::rotate_right_inplace(&mut aux_m_col[0], 1);
    utils::rotate_right_inplace(&mut aux_m_col[1], 2);
    utils::rotate_right_inplace(&mut aux_m_col[2], 3);
    utils::rotate_right_inplace(&mut aux_m_col[3], 3);

    for k in 0..4 {
        for round in 0..R - 2 {
            for i in 0..16 {
                let pos = 16 * round + i;
                let ys_pos = 16 * round + aux_m_col[k][i];
                let ys_offset = if k < 3 {
                    registry.s_box
                } else {
                    registry.m_col[0]
                };
                let v_even = v[(16 * (R - 2) * k + pos) * 2];
                let v_odd = v[(16 * (R - 2) * k + pos) * 2 + 1];
                dst[(registry.m_col[k] + pos) * 2] += v_even;
                dst[(ys_offset + ys_pos) * 2] += r * v_even;
                dst[(registry.m_col[k + 1] + pos) * 2] += r2 * v_even;

                dst[(registry.m_col[k] + pos) * 2 + 1] += v_odd;
                dst[(ys_offset + ys_pos) * 2 + 1] += r * v_odd;
                dst[(registry.m_col[k + 1] + pos) * 2 + 1] += r2 * v_odd;
            }
        }
    }
}

// pub fn vec_cipher_mcol<F: Field, const R: usize>(
//     c_xor1: F,
//     c_xor2: F,
// ) -> (Vec<F>, Vec<usize>, Vec<usize>) {
//     let identity = (0..16).collect::<Vec<_>>();
//     let regions = registry::aes_offsets::<R>();

//     let mut aux_m_col = vec![identity; 4];
//     crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[0], 1);
//     crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[1], 2);
//     crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[2], 3);
//     crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[3], 3);

//     let mut v: Vec<F> = Vec::new();
//     let mut idx: Vec<usize> = Vec::new();
//     let mut round_state: Vec<usize> = Vec::new();

//     let mut counter = 0;

//     for k in 0..4 {
//         for round in 0..R - 2 {
//             let row_offset = 16 * round;
//             let lhs_offset = regions.s_box;
//             let rhs_offset = regions.m_col[0];

//             let round_matrix = mcol_round_constrain(row_offset, lhs_offset, rhs_offset, output_offset, c);
// for i in 0..16 {
//     let pos = 16 * round + i;
//     let ys_pos = 16 * round + aux_m_col[k][i];
//     let ys_offset = if k < 3 {
//         regions.s_box
//     } else {
//         regions.m_col[0]
//     };

//     let input_l = (regions.m_col[k] + pos) * 2;
//     let input_r = (regions.m_col[k + 1] + pos) * 2;
//     let outout = (ys_offset + ys_pos) * 2;

//     v.push(F::ONE);
//     v.push(c_xor2);
//     v.push(F::ONE);
//     v.push(c_xor2);
//     v.push(c_xor1);
//     v.push(c_xor1);

//     idx.push(input_l);
//     idx.push(input_r);
//     idx.push(input_l + 1);
//     idx.push(input_r + 1);
//     idx.push(outout);
//     idx.push(outout + 1);

//     let c = [counter; 6];
//     round_state.extend_from_slice(&c);

//     counter += 1;
// }
//         }
//     }
//     (v, idx, round_state)
// }

pub fn gcm_final_xor<F: Field, const R: usize>(dst: &mut [F], v: &[F], r: F, r2: F) {
    let reg = registry::aes_gcm_block_offsets::<R>();
    let mut v_pos = 0;

    // XOR over 16 byte messages
    for i in 0..16 {
        let x_pos = 16 * i;
        let y_pos = 16 * (i + 1);
        let z_pos = 16 * (i + 2);

        let v_even = v[v_pos * 2];
        let v_odd = v[v_pos * 2 + 1];

        dst[(reg.final_xor + x_pos) * 2] += v_even;
        dst[(reg.final_xor + y_pos) * 2] += r * v_even;
        dst[(reg.final_xor + z_pos) * 2] += r2 * v_even;

        dst[(reg.final_xor + x_pos) * 2 + 1] += v_odd;
        dst[(reg.final_xor + y_pos) * 2 + 1] += r * v_odd;
        dst[(reg.final_xor + z_pos) * 2 + 1] += r2 * v_odd;

        v_pos += 1;
    }
}

pub fn cipher_addroundkey<F: Field, const R: usize>(
    output: &[u8; 16],
    dst: &mut [F],
    v: &[F],
    r: F,
    r2: F,
) -> F {
    let mut constant_term = F::from(0);
    let registry = registry::aes_offsets::<R>();

    for round in 0..R - 2 {
        for i in 0..16 {
            let pos = 16 * round + i;
            let v_even = v[pos * 2];
            let v_odd = v[pos * 2 + 1];
            dst[(registry.m_col[4] + pos) * 2] += v_even;
            dst[(registry.start + pos + 16) * 2] += r2 * v_even;
            dst[(registry.round_keys + pos + 16) * 2] += r * v_even;

            dst[(registry.m_col[4] + pos) * 2 + 1] += v_odd;
            dst[(registry.start + pos + 16) * 2 + 1] += r2 * v_odd;
            dst[(registry.round_keys + pos + 16) * 2 + 1] += r * v_odd;
        }
    }
    // final round
    #[allow(clippy::needless_range_loop)]
    for i in 0..16 {
        let pos = 16 * (R - 2) + i;
        let v_even = v[pos * 2];
        let v_odd = v[pos * 2 + 1];
        dst[(registry.s_box + pos) * 2] += v_even;
        dst[(registry.s_box + pos) * 2 + 1] += v_odd;
        dst[(registry.round_keys + pos + 16) * 2] += r * v_even;
        dst[(registry.round_keys + pos + 16) * 2 + 1] += r * v_odd;
        // in AES-EM mode, we would have to add the message instead.
        // dst[(OFFSETS.message + i) * 2] += r * v_even;
        // dst[(OFFSETS.message + i) * 2 + 1] += r * v_odd;
        constant_term += r2 * v_even * F::from(output[i] & 0xf);
        constant_term += r2 * v_odd * F::from(output[i] >> 4);
    }

    // initial round
    for i in 0..16 {
        let pos = 16 * (R - 1) + i;
        let v_even = v[pos * 2];
        let v_odd = v[pos * 2 + 1];
        // message
        dst[(registry.message + i) * 2] += v_even;
        dst[(registry.message + i) * 2 + 1] += v_odd;
        // initial round key
        dst[(registry.round_keys + i) * 2] += r * v_even;
        dst[(registry.round_keys + i) * 2 + 1] += r * v_odd;
        // .start
        dst[(registry.start + i) * 2] += r2 * v_even;
        dst[(registry.start + i) * 2 + 1] += r2 * v_odd;
    }
    constant_term
}

pub fn ks_lin_sbox_map<F: Field, const R: usize, const N: usize>(dst: &mut [F], v: &[F], r: F) {
    let reg = registry::aes_keysch_offsets::<R, N>();
    let n_4 = N / 4;
    let identity = [0, 1, 2, 3];
    let mut rotated_left = identity;
    rotated_left.rotate_left(1);

    for round in n_4..R {
        let idx = if N > 6 && (round * 4) % N == 4 {
            identity
        } else {
            rotated_left
        };
        for (y_j, x_j) in idx.into_iter().enumerate() {
            let x_pos = 16 * (round - n_4) + 3 * 4 + x_j;
            let y_pos = 4 * round + y_j;

            let c_lo = v[(round - n_4) * 4 + y_j];
            let c_hi = c_lo.double().double().double().double();
            dst[(reg.round_keys + x_pos) * 2] += c_lo;
            dst[(reg.round_keys + x_pos) * 2 + 1] += c_hi;
            dst[(reg.s_box + y_pos) * 2] += r * c_lo;
            dst[(reg.s_box + y_pos) * 2 + 1] += r * c_hi;
        }
    }
}

pub fn ks_lin_xor_map<F: Field, const R: usize, const N: usize>(
    dst: &mut [F],
    v: &[F],
    [r, r2]: [F; 2],
) -> F {
    let reg = registry::aes_keysch_offsets::<R, N>();
    // the running index over the source vector
    let mut v_pos = 0;
    // XXX. constant_term has to be mutated for supporting aes256 keyschedule
    let constant_term = F::from(0);

    // round_keys[i - n_4][1..4] XOR round_keys[i][0..3] = round_keys[i][1..4]
    let n_4 = N / 4;
    for round in n_4..R {
        for i in 1..4 {
            for j in 0..4 {
                let x_pos = 16 * (round - n_4) + i * 4 + j;
                let y_pos = 16 * round + (i - 1) * 4 + j;
                let z_pos = 16 * round + i * 4 + j;

                let v_even = v[v_pos * 2];
                let v_odd = v[v_pos * 2 + 1];

                dst[(reg.round_keys + x_pos) * 2] += v_even;
                dst[(reg.round_keys + y_pos) * 2] += r * v_even;
                dst[(reg.round_keys + z_pos) * 2] += r2 * v_even;

                dst[(reg.round_keys + x_pos) * 2 + 1] += v_odd;
                dst[(reg.round_keys + y_pos) * 2 + 1] += r * v_odd;
                dst[(reg.round_keys + z_pos) * 2 + 1] += r2 * v_odd;

                v_pos += 1;
            }
        }
    }

    // at this point,
    // v_pos = 3 * (R-1) * 4

    for round in n_4..R {
        for j in 0..4 {
            let x_pos = 16 * (round - n_4) + j;
            let y_pos = 4 * round + j;
            let z_pos = 16 * round + j;

            let v_even = v[v_pos * 2];
            let v_odd = v[v_pos * 2 + 1];

            dst[(reg.round_keys + x_pos) * 2] += v_even;
            dst[(reg.xor + y_pos) * 2] += r * v_even;
            dst[(reg.round_keys + z_pos) * 2] += r2 * v_even;

            dst[(reg.round_keys + x_pos) * 2 + 1] += v_odd;
            dst[(reg.xor + y_pos) * 2 + 1] += r * v_odd;
            dst[(reg.round_keys + z_pos) * 2 + 1] += r2 * v_odd;

            v_pos += 1;
        }
    }

    // at this point,
    // count = 3 * (R-1) * 4 + (R-1) * 4
    constant_term
}

pub fn sbox_constrain<F: Field, const R: usize>(c: F) -> linalg::SparseMatrix<F> {
    let reg = registry::aes_offsets::<R>();
    let sbox_mat = (0..R - 1)
        .map(|round| {
            let row_offset = round*16;
            let input_offset = reg.start + row_offset;
            let output_offset = reg.s_box + row_offset;
            sbox_round_constrain(row_offset, input_offset, output_offset, c)
        })
        .reduce(linalg::SparseMatrix::combine)
        .unwrap();
    sbox_mat
}

pub fn sbox_round_constrain<F: Field>(
    row_offset: usize,
    input_offset: usize,
    output_offset: usize,
    c: F,
) -> linalg::SparseMatrix<F> {
    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let s_row = utils::shiftrows(identity);
    let high = F::from(16);

    let rows = (0..16)
        .flat_map(|i| [row_offset + i; 4])
        .collect::<Vec<_>>();
    let cols = (0..16)
        .flat_map(|i| {
            [
                (s_row[i] + input_offset) * 2,
                (s_row[i] + input_offset) * 2 + 1,
                (i + output_offset) * 2,
                (i + output_offset) * 2 + 1,
            ]
        })
        .collect::<Vec<_>>();
    let vals = (0..16)
        .flat_map(|_| vec![F::ONE, high, c, c * high])
        .collect::<Vec<_>>();
    return linalg::SparseMatrix {
        num_rows: 16,
        vals,
        rows,
        cols,
    };
}

#[test]
fn test_sbox_constrain() {
    type F = ark_curve25519::Fr;
    use crate::traits::Witness;
    use ark_std::{UniformRand, Zero};
    use rand::Rng;

    let rng = &mut rand::thread_rng();
    let c = rng.gen();
    let mut sbox_mat: SparseMatrix<F> = sbox_constrain::<F, 11>(c);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let witness = crate::witness::cipher::AesCipherWitness::<F, 11, 4>::new(
        state,
        &key,
        F::zero(),
        F::zero(),
    );

    let vector_witness =
        crate::witness::cipher::AesCipherWitness::<F, 11, 4>::full_witness(&witness);
    let needles = &sbox_mat * &vector_witness;
    let haystack = haystack_sbox(c);

    assert!(
        needles.iter().all(|x| haystack.contains(x)),
        "Witness: {:?}, Needles: {:?}",
        &vector_witness,
        &needles
    );
}

#[test]
fn test_sbox_round_constrain() {
    use ark_curve25519::Fr as FF;
    use rand::Rng;

    let rng = &mut rand::thread_rng();
    let c = rng.gen();
    let mut sbox_mat: SparseMatrix<FF> = sbox_round_constrain(0, 0, 16, c);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let mut witness_u8: Vec<u8> = Vec::new();
    witness_u8.extend_from_slice(&state);
    witness_u8.extend_from_slice(&round_trace.s_box);
    witness_u8 = witness_u8.iter().flat_map(|x| [x & 0xf, x >> 4]).collect();

    let witness = witness_u8.iter().map(|x| FF::from(*x)).collect::<Vec<_>>();
    let needles = &sbox_mat * &witness;
    let haystack = haystack_sbox(c);

    assert!(
        needles.iter().all(|x| haystack.contains(x)),
        "Witness: {:?}, Needles: {:?}, Needles not in stack {:?}",
        &witness,
        &needles,
        &needles
            .iter()
            .filter(|&x| !haystack.contains(&x))
            .cloned()
            .collect::<Vec<_>>()
    );

    println!("Witness: {:?}, Needles: {:?}", &witness, &needles);
}

pub fn rj2_constrain<F: Field, const R: usize>(c: F) -> linalg::SparseMatrix<F> { 
    let reg = registry::aes_offsets::<R>();
    let rj2_mat = (0..R - 2)
        .map(|round| {
            let row_offset = round*16;
            let input_offset = reg.s_box + row_offset;
            let output_offset = reg.m_col[0] + row_offset;
            rj2_round_constrain(row_offset, input_offset, output_offset, c)
        })
        .reduce(linalg::SparseMatrix::combine)
        .unwrap();
    rj2_mat
}

fn rj2_round_constrain<F: Field>(
    row_offset: usize,
    input_offset: usize,
    output_offset: usize,
    c: F,
) -> linalg::SparseMatrix<F> {
    let high = F::from(16);
    let rows = (0..16)
        .flat_map(|i| [row_offset + i; 4])
        .collect::<Vec<_>>();
    let cols = (0..16)
        .map(|i| [i + input_offset, i + output_offset])
        .flat_map(|x| [x[0]*2, x[0]*2+1, x[1]*2, x[1]*2+1])
        .collect::<Vec<_>>();
    let vals = (0..16).flat_map(|_| vec![F::ONE, high, c, c*high]).collect::<Vec<_>>();
    return linalg::SparseMatrix {
        num_rows: 16,
        vals,
        rows,
        cols,
    };
}

#[test]
fn test_rj2_constrain() {
    type F = ark_curve25519::Fr;
    use crate::traits::Witness;
    use ark_std::{UniformRand, Zero};
    use rand::Rng;

    let rng = &mut rand::thread_rng();
    let c = rng.gen();
    let mut rj2_mat: SparseMatrix<F> = rj2_constrain::<F, 11>(c);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let witness = crate::witness::cipher::AesCipherWitness::<F, 11, 4>::new(
        state,
        &key,
        F::zero(),
        F::zero(),
    );

    let vector_witness =
        crate::witness::cipher::AesCipherWitness::<F, 11, 4>::full_witness(&witness);
    let needles = &rj2_mat * &vector_witness;
    let haystack = haystack_rj2(c);

    assert!(
        needles.iter().all(|x| haystack.contains(x)),
        "Witness: {:?}, Needles: {:?}",
        &vector_witness,
        &needles
    );
}

#[test]
fn test_rj2_round_constrain() {
    use ark_curve25519::Fr as FF;
    use rand::Rng;

    let rng = &mut rand::thread_rng();
    let c = rng.gen();
    let rj2_mat = rj2_round_constrain(0, 0, 16, c);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let mut witness_u8: Vec<u8> = Vec::new();
    witness_u8.extend_from_slice(&round_trace.s_box);
    witness_u8.extend_from_slice(&round_trace.m_col[0]);
    witness_u8 = witness_u8.iter().flat_map(|x| [x & 0xf, x >> 4]).collect();


    let witness = witness_u8.iter().map(|x| FF::from(*x)).collect::<Vec<_>>();
    let needles = &rj2_mat * &witness;
    let haystack = haystack_rj2(c);

    assert!(
        needles.iter().all(|x| haystack.contains(x)),
        "Witness: {:?}, Needles: {:?}",
        &witness,
        &needles
    );
}

fn mcol_round_constrain<F: Field>(
    row_offset: usize,
    lhs_offset: usize,
    rhs_offset: usize,
    output_offset: usize,
    c: F,
    c2: F
) -> linalg::SparseMatrix<F> {
    let mut m_col_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    rotate_right_inplace(&mut m_col_idx, 1);

    let rows = (0..32)
        .flat_map(|i| [row_offset + i; 3])
        .collect::<Vec<_>>();
    let cols = (0..16)
        .map(|i| [i + lhs_offset, m_col_idx[i] + rhs_offset, i + output_offset]).flat_map(|x| [x[0]*2, x[1]*2, x[2]*2,x[0]*2+1, x[1]*2+1, x[2]*2+1])
        .collect::<Vec<_>>();
    let vals = (0..32)
        .flat_map(|_| vec![F::ONE, c, c2])
        .collect::<Vec<_>>();
    return linalg::SparseMatrix {
        num_rows: 32,
        vals,
        rows,
        cols,
    };
}

#[test]
fn test_mcol_round_constrain() {
    use ark_curve25519::Fr as FF;
    use rand::Rng;

    let rng = &mut rand::thread_rng();
    let c = rng.gen::<FF>();
    let c2 = c.square();
    let mcol_mat = mcol_round_constrain(0, 16, 0, 32, c, c2);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let mut witness_u8: Vec<u8> = Vec::new();
    witness_u8.extend_from_slice(&round_trace.s_box);
    witness_u8.extend_from_slice(&round_trace.m_col[0]);
    witness_u8.extend_from_slice(&round_trace.m_col[1]);
    witness_u8 = witness_u8.iter().flat_map(|x| [x & 0xf, x >> 4]).collect();
    let witness = witness_u8.iter().map(|x| FF::from(*x)).collect::<Vec<_>>();


    // witness_u8.extend_from_slice(&round_trace.m_col[2]);
    // witness_u8.extend_from_slice(&round_trace.m_col[3]);

    // let witness_lo = witness_u8
    //     .iter()
    //     .map(|x| FF::from(*x & 0x0f))
    //     .collect::<Vec<_>>();
    // let witness_hi = witness_u8
    //     .iter()
    //     .map(|x| FF::from(*x >> 4))
    //     .collect::<Vec<_>>();

    // let needles_lo = &mcol_mat * &witness_lo;
    // let needles_hi = &mcol_mat * &witness_hi;
    // let needles = needles_lo
    //     .iter()
    //     .chain(needles_hi.iter())
    //     .collect::<Vec<_>>();
    let needles = &mcol_mat * &witness;
    let haystack = haystack_xor(c, c2);

    // assert!(needles.iter().all(|x| haystack.contains(x)),);

    assert!(
        needles.iter().all(|x| haystack.contains(x)),
        "Needles: {:?}, Needles not in stack {:?}",
        &needles,
        &needles
            .iter()
            .filter(|&x| !haystack.contains(&x))
            .cloned()
            .collect::<Vec<_>>()
    );
}

/// Return the lookup table for the AES s-box
#[inline]
fn haystack_sbox<F: Field>(c_sbox: F) -> Vec<F> {
    (0u8..=255)
        .map(|i| {
            let x = i;
            let y = utils::SBOX[x as usize];
            F::from(x) + c_sbox * F::from(y)
        })
        .collect()
}

#[inline]
fn haystack_rj2<F: Field>(c_rj2: F) -> Vec<F> {
    (0u8..=255)
        .map(|i| {
            let x = i;
            let y = utils::RJ2[x as usize];
            F::from(x) + c_rj2 * F::from(y)
        })
        .collect()
}

#[inline]
fn haystack_xor<F: Field>(c: F, c2: F) -> Vec<F> {
    (0u8..=0xff)
        .map(|i| {
            let x = i & 0x0f;
            let y = i >> 4;
            let z = x ^ y;
            F::from(x) + c * F::from(y) + c2 * F::from(z)
        })
        .collect()
}