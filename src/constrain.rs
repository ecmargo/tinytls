use crate::witness::{registry, trace::utils};
use ark_ff::Field;

pub fn aes_trace_to_needles<F: Field, const R: usize>(
    output: &[u8; 16],
    src: &[F],
    [c_xor, c_xor2, c_sbox, c_rj2]: [F; 4],
) -> (Vec<F>, F) {
    let reg = registry::aes_offsets::<R>();
    let mut dst = vec![F::ZERO; reg.witness_len * 2];
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
    let mut dst = vec![F::ZERO; registry.witness_len * 2];
    let mut offset: usize = 0;
    crate::constrain::ks_lin_sbox_map::<F, R, N>(&mut dst, src, c_sbox);
    offset += 4 * (R - N / 4);
    let constant_term =
        crate::constrain::ks_lin_xor_map::<F, R, N>(&mut dst, &src[offset..], [c_xor, c_xor2]);
    (dst, constant_term)
}

pub fn cipher_sbox<F: Field, const R: usize>(dst: &mut [F], v: &[F], r: F) {
    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let s_row = utils::shiftrows(identity);
    let reg = registry::aes_offsets::<R>();

    for round in 0..R - 1 {
        for i in 0..16 {
            let s_row_pos = 16 * round + s_row[i] as usize;
            let s_box_pos = 16 * round + i;
            let c_lo = v[round * 16 + i];
            let c_hi = c_lo.double().double().double().double();
            dst[(reg.start + s_row_pos) * 2] += c_lo;
            dst[(reg.start + s_row_pos) * 2 + 1] += c_hi;
            dst[(reg.s_box + s_box_pos) * 2] += r * c_lo;
            dst[(reg.s_box + s_box_pos) * 2 + 1] += r * c_hi;
        }
    }
}

pub fn cipher_rj2<F: Field, const R: usize>(dst: &mut [F], v: &[F], r: F) {
    let reg = registry::aes_offsets::<R>();

    for round in 0..R - 2 {
        for i in 0..16 {
            let pos = 16 * round + i;
            let c_lo = v[pos];
            let c_hi = c_lo.double().double().double().double();
            dst[(reg.s_box + pos) * 2] += c_lo;
            dst[(reg.s_box + pos) * 2 + 1] += c_hi;
            dst[(reg.m_col[0] + pos) * 2] += r * c_lo;
            dst[(reg.m_col[0] + pos) * 2 + 1] += r * c_hi;
        }
    }
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

// New functions that will need to be integrated above
#[cfg(test)]
pub fn xor_constrain_no_output<F: Field>(
    lhs_offset: usize,
    rhs_offset: usize,
    c: F,
) -> SparseMatrix<F> {
    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];

    let rows = (0..32).flat_map(|i| [i; 2]).collect::<Vec<_>>();
    let cols = (0..16)
        .map(|i| [i + lhs_offset, identity[i] as usize + rhs_offset])
        .flat_map(|x| [x[0] * 2, x[1] * 2, x[0] * 2 + 1, x[1] * 2 + 1])
        .collect::<Vec<_>>();
    let vals = (0..32).flat_map(|_| vec![F::ONE, c]).collect::<Vec<_>>();
    return SparseMatrix {
        num_rows: 32,
        vals,
        rows,
        cols,
    };
}

//Generic for generating constraints for XOR
#[cfg(test)]
pub fn rotate_xor_contstrain<F: Field>(
    lhs_offset: usize,
    rhs_offset: usize,
    output_offset: usize,
    c: F,
    c2: F,
    n_rotate: usize,
) -> SparseMatrix<F> {
    use crate::witness::trace::utils::rotate_right_inplace;

    let mut identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    rotate_right_inplace(&mut identity, n_rotate);

    let rows = (0..32).flat_map(|i| [i; 3]).collect::<Vec<_>>();
    let cols = (0..16)
        .map(|i| {
            [
                i + lhs_offset,
                identity[i] as usize + rhs_offset,
                i + output_offset,
            ]
        })
        .flat_map(|x| {
            [
                x[0] * 2,
                x[1] * 2,
                x[2] * 2,
                x[0] * 2 + 1,
                x[1] * 2 + 1,
                x[2] * 2 + 1,
            ]
        })
        .collect::<Vec<_>>();
    let vals = (0..32)
        .flat_map(|_| vec![F::ONE, c, c2])
        .collect::<Vec<_>>();
    return SparseMatrix {
        num_rows: 32,
        vals,
        rows,
        cols,
    };
}

//Generate lookup constraints for GCM XOR
#[cfg(test)]
pub fn gcm_final_xor_block_constrain<F: Field, const R: usize>(c: F) -> SparseMatrix<F> {
    let reg = registry::aes_gcm_block_offsets::<R>();
    let lhs_offset = reg.aes_output;
    let rhs_offset = reg.plain_text;

    let xor_mat = xor_constrain_no_output(lhs_offset, rhs_offset, c);

    xor_mat
}

//Generate constraints for a single add_roundkey
#[cfg(test)]
pub fn add_roundkey_round_constrain<F: Field>(
    lhs_offset: usize,
    rhs_offset: usize,
    output_offset: usize,
    c: F,
    c2: F,
) -> SparseMatrix<F> {
    rotate_xor_contstrain(lhs_offset, rhs_offset, output_offset, c, c2, 0)
}

//Generate all constraints for adding round key, specific for basic AES
#[cfg(test)]
pub fn add_roundkey_constrain_aes<F: Field, const R: usize>(c: F, c2: F) -> SparseMatrix<F> {
    //Normal rounds
    let reg = registry::aes_offsets::<R>();

    let mut add_roundkey_mat = SparseMatrix::new();

    let middle_rounds = (0..R - 2)
        .map(|round| {
            let offset_shift = round * 16;
            let lhs_offset = reg.m_col[4] + offset_shift;
            //Need an initial shift in rhs nand output offset because for round i we need the i+1 round key and start
            let rhs_offset = reg.round_keys + offset_shift + 16;
            let output_offset = reg.start + offset_shift + 16;

            add_roundkey_round_constrain::<F>(lhs_offset, rhs_offset, output_offset, c, c2)
        })
        .reduce(SparseMatrix::combine_with_rowshift)
        .unwrap();

    add_roundkey_mat = add_roundkey_mat.combine_with_rowshift(middle_rounds);

    let initial_round = {
        let lhs_offset = reg.message;
        let rhs_offset = reg.round_keys;
        let output_offset = reg.start;

        add_roundkey_round_constrain::<F>(lhs_offset, rhs_offset, output_offset, c, c2)
    };
    add_roundkey_mat = add_roundkey_mat.combine_with_rowshift(initial_round);

    let final_round = {
        let offset_shift = (R - 2) * 16;
        let lhs_offset = reg.s_box + offset_shift;
        let rhs_offset = reg.round_keys + offset_shift + 16;

        xor_constrain_no_output::<F>(lhs_offset, rhs_offset, c)
    };
    add_roundkey_mat = add_roundkey_mat.combine_with_rowshift(final_round);

    add_roundkey_mat
}

#[cfg(test)]
#[allow(unused)]
pub fn constant_term_constrain<F: Field>(c2: F) -> SparseMatrix<F> {
    let rows = (0..32).flat_map(|i| [i; 2]).collect::<Vec<_>>();
    let cols = (0..16)
        .map(|i| [i])
        .flat_map(|x| [x[0] * 2, x[0] * 2 + 1])
        .collect::<Vec<_>>();
    let vals = (0..32).flat_map(|_| vec![c2]).collect::<Vec<_>>();
    return SparseMatrix {
        num_rows: 32,
        vals,
        rows,
        cols,
    };
}

//Generate constraints for a single sbox round
#[cfg(test)]
pub fn sbox_round_constrain<F: Field>(
    input_offset: usize,
    output_offset: usize,
    c: F,
) -> SparseMatrix<F> {
    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let s_row = utils::shiftrows(identity);
    let high = F::from(16);

    let rows = (0..16).flat_map(|i| [i; 4]).collect::<Vec<_>>();
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
    return SparseMatrix {
        num_rows: 16,
        vals,
        rows,
        cols,
    };
}

//Generate constraints for all sbox rounds
#[cfg(test)]
pub fn sbox_constrain<F: Field, const R: usize>(c: F) -> SparseMatrix<F> {
    let reg = registry::aes_offsets::<R>();
    let sbox_mat = (0..R - 1)
        .map(|_| {
            let input_offset = reg.start;
            let output_offset = reg.s_box;
            sbox_round_constrain(input_offset, output_offset, c)
        })
        .reduce(SparseMatrix::combine_with_rowshift)
        .unwrap();
    sbox_mat
}

//Generate constraints for a single rj2 round
#[cfg(test)]
fn rj2_round_constrain<F: Field>(
    input_offset: usize,
    output_offset: usize,
    c: F,
) -> SparseMatrix<F> {
    let high = F::from(16);
    let rows = (0..16).flat_map(|i| [i; 4]).collect::<Vec<_>>();
    let cols = (0..16)
        .map(|i| [i + input_offset, i + output_offset])
        .flat_map(|x| [x[0] * 2, x[0] * 2 + 1, x[1] * 2, x[1] * 2 + 1])
        .collect::<Vec<_>>();
    let vals = (0..16)
        .flat_map(|_| vec![F::ONE, high, c, c * high])
        .collect::<Vec<_>>();
    return SparseMatrix {
        num_rows: 16,
        vals,
        rows,
        cols,
    };
}

//Generate constraints for all rj2 rounds
#[cfg(test)]
pub fn rj2_constrain<F: Field, const R: usize>(c: F) -> SparseMatrix<F> {
    let reg = registry::aes_offsets::<R>();
    let rj2_mat = (0..R - 2)
        .map(|_| {
            let input_offset = reg.s_box;
            let output_offset = reg.m_col[0];
            rj2_round_constrain(input_offset, output_offset, c)
        })
        .reduce(SparseMatrix::combine_with_rowshift)
        .unwrap();
    rj2_mat
}

//Generate constraints for all rounds of mcol
#[cfg(test)]
pub fn mcol_constrain<F: Field, const R: usize>(c: F, c2: F) -> SparseMatrix<F> {
    let reg = registry::aes_offsets::<R>();

    let mut mcol_mat = SparseMatrix::new();

    let mcol_mat1s: SparseMatrix<F> = (0..R - 2)
        .map(|round| {
            let offset_shift = 16 * round;
            let lhs_offset = reg.m_col[0] + offset_shift;
            let rhs_offset = reg.s_box + offset_shift;
            let output_offset = reg.m_col[1] + offset_shift;
            rotate_xor_contstrain(lhs_offset, rhs_offset, output_offset, c, c2, 1)
        })
        .reduce(SparseMatrix::combine_with_rowshift)
        .unwrap();

    mcol_mat = mcol_mat.combine(mcol_mat1s);

    let mcol_mat2s: SparseMatrix<F> = (0..R - 2)
        .map(|round| {
            let offset_shift = 16 * round;
            let lhs_offset = reg.m_col[1] + offset_shift;
            let rhs_offset = reg.s_box + offset_shift;
            let output_offset = reg.m_col[2] + offset_shift;
            rotate_xor_contstrain(lhs_offset, rhs_offset, output_offset, c, c2, 2)
        })
        .reduce(SparseMatrix::combine_with_rowshift)
        .unwrap();

    mcol_mat = mcol_mat.combine_with_rowshift(mcol_mat2s);

    let mcol_mat3s: SparseMatrix<F> = (0..R - 2)
        .map(|round| {
            let offset_shift = 16 * round;
            let lhs_offset = reg.m_col[2] + offset_shift;
            let rhs_offset = reg.s_box + offset_shift;
            let output_offset = reg.m_col[3] + offset_shift;
            rotate_xor_contstrain(lhs_offset, rhs_offset, output_offset, c, c2, 3)
        })
        .reduce(SparseMatrix::combine_with_rowshift)
        .unwrap();

    mcol_mat = mcol_mat.combine_with_rowshift(mcol_mat3s);

    let mcol_mat4s: SparseMatrix<F> = (0..R - 2)
        .map(|round| {
            let offset_shift = 16 * round;
            let lhs_offset = reg.m_col[3] + offset_shift;
            let rhs_offset = reg.m_col[0] + offset_shift;
            let output_offset = reg.m_col[4] + offset_shift;
            rotate_xor_contstrain(lhs_offset, rhs_offset, output_offset, c, c2, 3)
        })
        .reduce(SparseMatrix::combine_with_rowshift)
        .unwrap();

    mcol_mat = mcol_mat.combine_with_rowshift(mcol_mat4s);

    mcol_mat
}

#[cfg(test)]
mod tests {
    use ark_ff::{AdditiveGroup, Field};
    use rand::Rng;

    use crate::linalg::SparseMatrix;
    use crate::witness::trace::utils;
    use crate::Witness;

    use super::*;

    /// Return the lookup table for the AES s-box
    fn haystack_sbox<F: Field>(c_sbox: F) -> Vec<F> {
        (0u8..=255)
            .map(|i| {
                let x = i;
                let y = utils::SBOX[x as usize];
                F::from(x) + c_sbox * F::from(y)
            })
            .collect()
    }

    /// Return the lookup table for the AES rj2
    fn haystack_rj2<F: Field>(c_rj2: F) -> Vec<F> {
        (0u8..=255)
            .map(|i| {
                let x = i;
                let y = utils::RJ2[x as usize];
                F::from(x) + c_rj2 * F::from(y)
            })
            .collect()
    }

    /// Return the lookup table for the AES XOR
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

    #[test]
    fn test_add_roundkey_round_constrain() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c: F = rng.gen();
        let c2 = c.square();

        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();

        let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
        let mut witness_u8: Vec<u8> = Vec::new();
        witness_u8.extend_from_slice(&round_trace.start);
        witness_u8.extend_from_slice(&round_trace.m_col[4]);
        witness_u8.extend_from_slice(&key);
        witness_u8 = witness_u8.iter().flat_map(|x| [x & 0xf, x >> 4]).collect();

        let witness = witness_u8.iter().map(|x| F::from(*x)).collect::<Vec<_>>();

        let add_roundkey_mat = add_roundkey_round_constrain(16, 32, 0, c, c2);

        let needles = &add_roundkey_mat * &witness;
        let haystack = haystack_xor(c, c2);

        assert!(
            needles.iter().all(|x| haystack.contains(x)),
            "Needles: {:?}, Needles not in stack {:?}",
            &needles.len(),
            &needles
                .iter()
                .filter(|&x| !haystack.contains(&x))
                .cloned()
                .collect::<Vec<_>>()
                .len()
        );
    }

    #[test]
    fn test_final_xor_constrain_gcm() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c: F = rng.gen();
        let c2 = c.square();

        let final_xor_mat: SparseMatrix<F> = gcm_final_xor_block_constrain::<F, 15>(c);
        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();
        let iv: [u8; 12] = rng.gen::<[u8; 12]>();
        let mut ctr = crate::witness::trace::gcm::AesGCMCounter::create_icb(iv);
        ctr.count += 1;

        let witness = crate::witness::gcm::AesGCMCipherBlockWitness::<F, 15, 8>::new(
            ctr,
            &key,
            state,
            F::ZERO,
            F::ZERO,
        );

        let vector_witness =
            crate::witness::gcm::AesGCMCipherBlockWitness::<F, 15, 8>::full_witness(&witness);
        let needles = &final_xor_mat * &vector_witness;
        let haystack = haystack_xor(c, c2);

        assert!(
            needles.iter().all(|x| haystack.contains(x)),
            "Needles: {:?}, Needles not in stack {:?}",
            &needles.len(),
            &needles
                .iter()
                .filter(|&x| !haystack.contains(&x))
                .cloned()
                .collect::<Vec<_>>()
                .len()
        );
    }

    #[test]
    fn test_add_roundkey_constrain_aes() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c: F = rng.gen();
        let c2 = c.square();

        let add_roundkey_mat: SparseMatrix<F> = add_roundkey_constrain_aes::<F, 11>(c, c2);
        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();
        let witness = crate::witness::cipher::AesCipherWitness::<F, 11, 4>::new(
            state,
            &key,
            F::ZERO,
            F::ZERO,
        );

        let vector_witness =
            crate::witness::cipher::AesCipherWitness::<F, 11, 4>::full_witness(&witness);
        let needles = &add_roundkey_mat * &vector_witness;
        let haystack = haystack_xor(c, c2);

        assert!(
            needles.iter().all(|x| haystack.contains(x)),
            "Needles: {:?}, Needles not in stack {:?}",
            &needles.len(),
            &needles
                .iter()
                .filter(|&x| !haystack.contains(&x))
                .cloned()
                .collect::<Vec<_>>()
                .len()
        );
    }

    #[test]
    fn test_sbox_round_constrain() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c = rng.gen();

        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();

        let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
        let mut witness_u8: Vec<u8> = Vec::new();
        witness_u8.extend_from_slice(&state);
        witness_u8.extend_from_slice(&round_trace.s_box);
        witness_u8 = witness_u8.iter().flat_map(|x| [x & 0xf, x >> 4]).collect();

        let witness = witness_u8.iter().map(|x| F::from(*x)).collect::<Vec<_>>();

        //Input offeset - state (start of witness_u8)
        //Output offset - round_trace.s_box (state + 16)
        let sbox_mat: SparseMatrix<F> = sbox_round_constrain(0, 16, c);

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
    }

    #[test]
    fn test_sbox_constrain() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c = rng.gen();

        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();

        let witness = crate::witness::cipher::AesCipherWitness::<F, 11, 4>::new(
            state,
            &key,
            F::ZERO,
            F::ZERO,
        );

        let vector_witness =
            crate::witness::cipher::AesCipherWitness::<F, 11, 4>::full_witness(&witness);

        let sbox_mat: SparseMatrix<F> = sbox_constrain::<F, 11>(c);

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
    #[ignore = "Not ready yet"]
    fn test_constant_constrain() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();

        let witness = crate::witness::cipher::AesCipherWitness::<F, 11, 4>::new(
            state,
            &key,
            F::ZERO,
            F::ZERO,
        );

        let vector_witness =
            crate::witness::cipher::AesCipherWitness::<F, 11, 4>::full_witness(&witness);

        // let output = witness
        //     .trace
        //     .output
        //     .iter()
        //     .flat_map(|x| [x & 0xf, x >> 4])
        //     .collect::<Vec<_>>();

        let c = rng.gen::<F>();
        // let c2 = rng.gen::<F>();
        // let constant_mat: SparseMatrix<F> = constant_term_constrain::<F>(c2);

        // let needles = constant_mat * &output;
        let needles = Vec::new();
        let haystack = haystack_sbox(c);
        assert!(
            needles.iter().all(|x| haystack.contains(x)),
            "Witness: {:?}, Needles: {:?}",
            &vector_witness,
            &needles
        );
    }

    #[test]
    fn test_rj2_round_constrain() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c = rng.gen();

        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();

        let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
        let mut witness_u8: Vec<u8> = Vec::new();
        witness_u8.extend_from_slice(&round_trace.s_box);
        witness_u8.extend_from_slice(&round_trace.m_col[0]);
        witness_u8 = witness_u8.iter().flat_map(|x| [x & 0xf, x >> 4]).collect();
        let witness = witness_u8.iter().map(|x| F::from(*x)).collect::<Vec<_>>();

        //Input offset - round_trace.s_box (start of witness_u8)
        //Ouput offset - round_trace.m_col[0](input + 16)
        let rj2_mat = rj2_round_constrain(0, 16, c);

        let needles = &rj2_mat * &witness;
        let haystack = haystack_rj2(c);

        assert!(
            needles.iter().all(|x| haystack.contains(x)),
            "Witness: {:?}, Needles: {:?}",
            &witness,
            &needles
        );
    }

    #[test]
    fn test_rj2_constrain() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c = rng.gen();

        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();

        let witness = crate::witness::cipher::AesCipherWitness::<F, 11, 4>::new(
            state,
            &key,
            F::ZERO,
            F::ZERO,
        );

        let vector_witness =
            crate::witness::cipher::AesCipherWitness::<F, 11, 4>::full_witness(&witness);

        let rj2_mat: SparseMatrix<F> = rj2_constrain::<F, 11>(c);
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
    fn test_mcol_constrain() {
        type F = ark_curve25519::Fr;

        let rng = &mut rand::thread_rng();
        let c = rng.gen::<F>();
        let c2 = c.square();

        let mcol_mat = mcol_constrain::<F, 11>(c, c2);

        let state = rng.gen::<[u8; 16]>();
        let key = rng.gen::<[u8; 16]>();

        let witness = crate::witness::cipher::AesCipherWitness::<F, 11, 4>::new(
            state,
            &key,
            F::ZERO,
            F::ZERO,
        );

        let vector_witness =
            crate::witness::cipher::AesCipherWitness::<F, 11, 4>::full_witness(&witness);
        let needles = &mcol_mat * &vector_witness;
        let haystack = haystack_xor(c, c2);

        assert!(
            needles.iter().all(|x| haystack.contains(x)),
            "Needles: {:?}, Needles not in stack {:?}",
            &needles.len(),
            &needles
                .iter()
                .filter(|&x| !haystack.contains(&x))
                .cloned()
                .collect::<Vec<_>>()
                .len()
        );
    }
}
