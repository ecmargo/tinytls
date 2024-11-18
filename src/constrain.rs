use crate::linalg;
use crate::witness::trace::utils::rotate_right_inplace;
use crate::witness::{registry, trace::utils};
use ark_ff::Field;

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

//Make sparse matrix struct (I feel like this must exist already in rust and adapt sbox to work for a single round on that)

pub fn combine_yale_to_needles<F: Field>(
    round_states: Vec<usize>,
    idx: Vec<usize>,
    selectors: Vec<F>,
    witness: Vec<F>,
) -> Vec<F> {
    let mut output: Vec<F> = Vec::new();

    for count in 0..round_states[round_states.len() - 1] {
        let round_indices: Vec<usize> = round_states
            .iter()
            .enumerate()
            .filter(|(_, &rs)| rs == count)
            .map(|(index, _)| index)
            .collect();

        let idxs: Vec<usize> = round_indices.iter().map(|&i| idx[i]).collect();

        let select: Vec<F> = round_indices.iter().map(|&i| selectors[i]).collect();
        let wits: Vec<F> = idxs.iter().map(|&i| witness[i]).collect();

        output.push(linalg::inner_product(&select, &wits));
    }
    output
}

pub fn vec_cipher_sbox<F: Field, const R: usize>(c_sbox: F) -> (Vec<F>, Vec<usize>, Vec<usize>) {
    let regions = registry::aes_offsets::<R>();

    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let s_row = utils::shiftrows(identity);

    let high = F::from(16);
    let mut input;
    let mut output;

    let mut v: Vec<F> = Vec::new();
    let mut cols: Vec<usize> = Vec::new();
    let mut rows: Vec<usize> = Vec::new();

    let mut counter = 0;

    //Structure as a Vec<Vec<F>> where you push after each inner iteration
    for round in 0..R - 1 {
        for i in 0..16 {
            let s_row_pos = 16 * round + s_row[i] as usize;
            let s_box_pos = 16 * round + i;
            input = (regions.start + s_row_pos) * 2;
            output = (regions.s_box + s_box_pos) * 2;

            v.push(F::ONE);
            v.push(high);
            v.push(c_sbox);
            v.push(c_sbox * high);

            cols.push(input);
            cols.push(input + 1);
            cols.push(output);
            cols.push(output + 1);

            let row: [usize; 4] = [counter; 4];
            rows.extend_from_slice(&row);

            counter += 1;
        }
    }
    (v, cols, rows)
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

pub fn vec_cipher_rj2<F: Field, const R: usize>(c_rj2: F) -> (Vec<F>, Vec<usize>, Vec<usize>) {
    let regions = registry::aes_offsets::<R>();

    let high = F::from(16);
    let mut input;
    let mut output;

    let mut v: Vec<F> = Vec::new();
    let mut idx: Vec<usize> = Vec::new();
    let mut round_state: Vec<usize> = Vec::new();

    let mut counter = 0;

    for round in 0..R - 2 {
        for i in 0..16 {
            let pos = 16 * round + i;
            input = (regions.s_box + pos) * 2;
            output = (regions.m_col[0] + pos) * 2;

            v.push(F::ONE);
            v.push(high);
            v.push(c_rj2);
            v.push(c_rj2 * high);

            idx.push(input);
            idx.push(input + 1);
            idx.push(output);
            idx.push(output + 1);

            let c = [counter; 4];
            round_state.extend_from_slice(&c);

            counter += 1;
        }
    }
    (v, idx, round_state)
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

pub fn vec_cipher_mcol<F: Field, const R: usize>(
    c_xor1: F,
    c_xor2: F,
) -> (Vec<F>, Vec<usize>, Vec<usize>) {
    let identity = (0..16).collect::<Vec<_>>();
    let regions = registry::aes_offsets::<R>();

    let mut aux_m_col = vec![identity; 4];
    crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[0], 1);
    crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[1], 2);
    crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[2], 3);
    crate::witness::trace::utils::rotate_right_inplace(&mut aux_m_col[3], 3);

    let mut v: Vec<F> = Vec::new();
    let mut idx: Vec<usize> = Vec::new();
    let mut round_state: Vec<usize> = Vec::new();

    let mut counter = 0;

    for k in 0..4 {
        for round in 0..R - 2 {
            for i in 0..16 {
                let pos = 16 * round + i;
                let ys_pos = 16 * round + aux_m_col[k][i];
                let ys_offset = if k < 3 {
                    regions.s_box
                } else {
                    regions.m_col[0]
                };

                let input_l = (regions.m_col[k] + pos) * 2;
                let input_r = (regions.m_col[k + 1] + pos) * 2;
                let outout = (ys_offset + ys_pos) * 2;

                v.push(F::ONE);
                v.push(c_xor2);
                v.push(F::ONE);
                v.push(c_xor2);
                v.push(c_xor1);
                v.push(c_xor1);

                idx.push(input_l);
                idx.push(input_r);
                idx.push(input_l + 1);
                idx.push(input_r + 1);
                idx.push(outout);
                idx.push(outout + 1);

                let c = [counter; 6];
                round_state.extend_from_slice(&c);

                counter += 1;
            }
        }
    }
    (v, idx, round_state)
}

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

fn sbox_round_constrain<F: Field>(
    row_offset: usize,
    input_offset: usize,
    output_offset: usize,
    c: F,
) -> linalg::SparseMatrix<F> {
    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let s_row = utils::shiftrows(identity);

    let rows = (0..16)
        .flat_map(|i| [row_offset + i; 2])
        .collect::<Vec<_>>();
    let cols = (0..16)
        .flat_map(|i| [s_row[i] + input_offset, i + output_offset])
        .collect::<Vec<_>>();
    let vals = (0..16).flat_map(|_| vec![F::ONE, c]).collect::<Vec<_>>();
    return linalg::SparseMatrix {
        num_rows: 16,
        vals,
        rows,
        cols,
    };
}

#[test]
fn test_sbox_constrain() {
    use ark_curve25519::Fr as FF;
    use rand::Rng;

    let rng = &mut rand::thread_rng();
    let c = rng.gen();
    let sbox_mat = sbox_round_constrain(0, 0, 16, c);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let mut witness_u8: Vec<u8> = Vec::new();
    witness_u8.extend_from_slice(&state);
    witness_u8.extend_from_slice(&round_trace.s_box);

    let witness = witness_u8.iter().map(|x| FF::from(*x)).collect::<Vec<_>>();
    let needles = &sbox_mat * &witness;
    let haystack = haystack_sbox(c);

    assert!(
        needles.iter().all(|x| haystack.contains(x)),
        "Witness: {:?}, Needles: {:?}",
        &witness,
        &needles
    );
}

fn rj2_round_constrain<F: Field>(
    row_offset: usize,
    input_offset: usize,
    output_offset: usize,
    c: F,
) -> linalg::SparseMatrix<F> {
    let rows = (0..16)
        .flat_map(|i| [row_offset + i; 2])
        .collect::<Vec<_>>();
    let cols = (0..16)
        .flat_map(|i| [i + input_offset, i + output_offset])
        .collect::<Vec<_>>();
    let vals = (0..16).flat_map(|_| vec![F::ONE, c]).collect::<Vec<_>>();
    return linalg::SparseMatrix {
        num_rows: 16,
        vals,
        rows,
        cols,
    };
}

#[test]
fn test_rj2_round_constrain() {
    use ark_curve25519::Fr as FF;
    use rand::Rng;

    let rng = &mut rand::thread_rng();
    let c = rng.gen();
    let sbox_mat = rj2_round_constrain(0, 0, 16, c);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let mut witness_u8: Vec<u8> = Vec::new();
    witness_u8.extend_from_slice(&round_trace.s_box);
    witness_u8.extend_from_slice(&round_trace.m_col[0]);

    let witness = witness_u8.iter().map(|x| FF::from(*x)).collect::<Vec<_>>();
    let needles = &sbox_mat * &witness;
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
) -> linalg::SparseMatrix<F> {
    let c2 = c.square();
    let mut m_col_idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    rotate_right_inplace(&mut m_col_idx, 1);

    let rows = (0..16)
        .flat_map(|i| [row_offset + i; 3])
        .collect::<Vec<_>>();
    let cols = (0..16)
        .flat_map(|i| [i + lhs_offset, m_col_idx[i] + rhs_offset, i + output_offset])
        .collect::<Vec<_>>();
    let vals = (0..16)
        .flat_map(|_| vec![F::ONE, c, c2])
        .collect::<Vec<_>>();
    return linalg::SparseMatrix {
        num_rows: 16,
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
    let mcol_mat = mcol_round_constrain(0, 16, 0, 32, c);
    let state = rng.gen::<[u8; 16]>();
    let key = rng.gen::<[u8; 16]>();
    let round_trace = crate::witness::trace::cipher::aes_round_trace(state, key);
    let mut witness_u8: Vec<u8> = Vec::new();
    witness_u8.extend_from_slice(&round_trace.s_box);
    witness_u8.extend_from_slice(&round_trace.m_col[0]);
    witness_u8.extend_from_slice(&round_trace.m_col[1]);
    // witness_u8.extend_from_slice(&round_trace.m_col[2]);
    // witness_u8.extend_from_slice(&round_trace.m_col[3]);


    let witness_lo = witness_u8.iter().map(|x| FF::from(*x & 0x0f)).collect::<Vec<_>>();
    let witness_hi = witness_u8.iter().map(|x| FF::from(*x >> 4)).collect::<Vec<_>>();

    let needles_lo = &mcol_mat * &witness_lo;
    let needles_hi = &mcol_mat * &witness_hi;
    let needles = needles_lo.iter().chain(needles_hi.iter()).collect::<Vec<_>>();
    let haystack = haystack_xor(c, c2);

    assert!(
        needles.iter().all(|x| haystack.contains(x)),
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

#[test]
fn test_sbox() {
    type F = ark_curve25519::Fr;
    use crate::traits::Witness;
    use crate::witness::cipher::AesCipherWitness;
    use ark_std::{UniformRand, Zero};

    let rng = &mut rand::thread_rng();

    let message = [
        0x4A, 0x8F, 0x6D, 0xE2, 0x12, 0x7B, 0xC9, 0x34, 0xA5, 0x58, 0x91, 0xFD, 0x23, 0x69, 0x0C,
        0xE7,
    ];
    let key = [
        0xE7u8, 0x4A, 0x8F, 0x6D, 0xE2, 0x12, 0x7B, 0xC9, 0x34, 0xA5, 0x58, 0x91, 0xFD, 0x23, 0x69,
        0x0C,
    ];

    let c_sbox = F::rand(rng);

    let witness = AesCipherWitness::<F, 11, 4>::new(message, &key, F::zero(), F::zero());

    let vector_witness = AesCipherWitness::<F, 11, 4>::full_witness(&witness);

    let (v, idx, round_state) = vec_cipher_sbox::<F, 11>(c_sbox);

    let output: Vec<F> = combine_yale_to_needles::<F>(round_state, idx, v, vector_witness);

    let haystack_s_box = haystack_sbox(c_sbox);

    assert!(output.into_iter().all(|x| (haystack_s_box.contains(&x))));
}

// #[test]
// fn test_rj2() {
//     type F = ark_curve25519::Fr;
//     use crate::traits::Witness;
//     use crate::witness::cipher::AesCipherWitness;
//     use ark_std::{UniformRand, Zero};

//     let rng = &mut rand::thread_rng();

//     let message = [
//         0x4A, 0x8F, 0x6D, 0xE2, 0x12, 0x7B, 0xC9, 0x34, 0xA5, 0x58, 0x91, 0xFD, 0x23, 0x69, 0x0C,
//         0xE7,
//     ];
//     let key = [
//         0xE7u8, 0x4A, 0x8F, 0x6D, 0xE2, 0x12, 0x7B, 0xC9, 0x34, 0xA5, 0x58, 0x91, 0xFD, 0x23, 0x69,
//         0x0C,
//     ];

//     let c_rj2 = F::rand(rng);

//     let haystack_r2j = (0u8..=255)
//         .map(|i| {
//             let x = i;
//             let y = utils::RJ2[x as usize];
//             F::from(x) + c_rj2 * F::from(y)
//         })
//         .collect::<Vec<_>>();

//     let witness = AesCipherWitness::<F, 11, 4>::new(message, &key, F::zero(), F::zero());

//     let vector_witness = AesCipherWitness::<F, 11, 4>::full_witness(&witness);

//     //Leaving these here because not entirely sure what to do with them.
//     let offset = 16 * 10;
//     let sub = (&vector_witness[offset..]).to_vec();

//     let (v, idx, round_state) = vec_cipher_rj2::<F, 11>(c_rj2);

//     let output: Vec<F> = combine_yale_to_needles(round_state, idx, v, vector_witness);

//     assert!(output.into_iter().all(|x| (haystack_r2j.contains(&x))));
// }

// #[test]
// fn test_mcol() {
//     use ark_ff::{AdditiveGroup, Field};
//     type F = ark_curve25519::Fr;
//     use crate::traits::Witness;
//     use crate::witness::cipher::AesCipherWitness;
//     use ark_std::{UniformRand, Zero};

//     let rng = &mut rand::thread_rng();

//     let message = [
//         0x4A, 0x8F, 0x6D, 0xE2, 0x12, 0x7B, 0xC9, 0x34, 0xA5, 0x58, 0x91, 0xFD, 0x23, 0x69, 0x0C,
//         0xE7,
//     ];
//     let key = [
//         0xE7u8, 0x4A, 0x8F, 0x6D, 0xE2, 0x12, 0x7B, 0xC9, 0x34, 0xA5, 0x58, 0x91, 0xFD, 0x23, 0x69,
//         0x0C,
//     ];

    // let c_xor1 = F::ZERO;
    // //F::rand(rng);
    // let c_xor2 = F::ZERO;
    // //F::rand(rng);


    // let haystack_xor = (0u8..=255)
    // .map(|i| {
    //     let x = i & 0xf;
    //     let y = i >> 4;
    //     let z = x ^ y;
    //     F::from(x) + c_xor1 * F::from(y) + c_xor2 * F::from(z)
    // })
    // .collect::<Vec<_>>();

    // let (haystack, inv_haystack) = lookup::compute_haystack([c_xor1, c_xor2, c_sbox, c_rj2], F::ZERO);

//     let witness = AesCipherWitness::<F, 11, 4>::new(message, &key, F::ZERO, F::ZERO);

//     let vector_witness = AesCipherWitness::<F, 11, 4>::full_witness(&witness);

    //Leaving these here because not entirely sure what to do with them.
    // let offset = 16 * 10;
    // let sub = (&vector_witness[offset..]).to_vec();

//     let (v, idx, round_state) = vec_cipher_mcol::<F, 11>(c_xor1, c_xor2);

//     let output: Vec<F> = combine_yale_to_needles(round_state, idx, v, vector_witness);

//     println!("output: {:?}", output); 
//     println!("haystack {:?}", haystack_xor);

//     assert!(output.into_iter().all(|x|(haystack_xor.contains(&x))));
// }

