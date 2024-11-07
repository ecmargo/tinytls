use ark_ff::{AdditiveGroup, Field};

use crate::{aes_utils, registry, Witness};

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

pub fn aes_gcm_block_trace_to_needles<F: Field, const R: usize>(
    aes_output: &[u8; 16],
    final_xor: &[u8; 16],
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

pub fn vec_cipher_sbox<F: Field, const R: usize>(c_sbox: F) -> (Vec<F>, Vec<usize>, Vec<usize>, usize) {
    let regions = registry::aes_offsets::<R>();

    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let s_row = aes_utils::shiftrows(identity);

    let high = F::from(16);
    let mut input;
    let mut output;

    let mut v: Vec<F> = Vec::new();
    let mut idx: Vec<usize> = Vec::new();
    let mut round_state: Vec<usize> = Vec::new();

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
            v.push(c_sbox*high);

            idx.push(input);
            idx.push(input + 1);
            idx.push(output);
            idx.push(output + 1);

            let c = [counter; 4];
            round_state.extend_from_slice(&c);

            counter += 1;
        }
    }
    (v, idx, round_state, counter)
}

#[test]
fn test_xi_sbox() {
    use crate::linalg;
    type F = ark_curve25519::Fr;
    use crate::witness_plain::AesCipherWitness;
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

    let wit_f: Vec<F> = vector_witness.into_iter().map(F::from).collect();

    let (v, idx, round_state, counter) = vec_cipher_sbox::<F, 11>(c_sbox);

    let mut output: Vec<F> = Vec::new();

    for count in 0..counter { 
        let round_indices: Vec<usize> = round_state.iter().enumerate().filter(|(_, &rs)| rs == count).map(|(index, _)| index).collect(); 

        let idxs: Vec<usize> = round_indices.iter().map(|&i| idx[i]).collect();

        let select: Vec<F> = round_indices.iter().map(|&i| v[i]).collect();
        let wits: Vec<F> = idxs.iter().map(|&i| wit_f[i]).collect();

        output.push(linalg::inner_product(&select, &wits));
    }

    let haystack_s_box = (0u8..=255)
        .map(|i| {
            let x = i;
            let y = aes_utils::SBOX[x as usize];
            F::from(x) + c_sbox * F::from(y)
        })
        .collect::<Vec<_>>();

    assert!(output.into_iter().all(|x|(haystack_s_box.contains(&x))));
}

pub fn cipher_sbox<F: Field, const R: usize>(dst: &mut [F], v: &[F], r: F) {
    let identity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15];
    let s_row = aes_utils::shiftrows(identity);
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

pub fn vec_cipher_rj2<F: Field, const R: usize>(c_rj2: F) -> (Vec<F>, Vec<usize>, Vec<usize>, usize) {
    let regions = registry::aes_offsets::<R>();

    let high = F::from(16);
    let mut input;
    let mut output;

    let mut v: Vec<F> = Vec::new();
    let mut idx: Vec<usize> = Vec::new();
    let mut round_state: Vec<usize> = Vec::new();

    let mut counter = 0;

    for round in 0..R - 1 {
        for i in 0..16 {
            let pos = 16 * round + i;
            input = (regions.s_box + pos) * 2;
            output = (regions.m_col[0] + pos) * 2;

            v.push(F::ONE); 
            v.push(high); 
            v.push(c_rj2);
            v.push(c_rj2*high);

            idx.push(input);
            idx.push(input + 1);
            idx.push(output);
            idx.push(output + 1);

            let c = [counter; 4];
            round_state.extend_from_slice(&c);

            counter += 1;
        }
    }
    (v, idx, round_state, counter)
}

#[test]
fn test_rj2() {
    use crate::linalg;
    type F = ark_curve25519::Fr;
    use crate::witness_plain::AesCipherWitness;
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

    let c_rj2 = F::ONE;
    //F::rand(rng);

    let haystack_r2j = (0u8..=255)
    .map(|i| {
        let x = i;
        let y = aes_utils::RJ2[x as usize];
        F::from(x) + c_rj2 * F::from(y)
    })
    .collect::<Vec<_>>();
    // println!("{:?}", haystack_r2j);

    let witness = AesCipherWitness::<F, 11, 4>::new(message, &key, F::zero(), F::zero());

    let vector_witness = AesCipherWitness::<F, 11, 4>::full_witness(&witness);

    let wit_f: Vec<F> = vector_witness.into_iter().map(F::from).collect();

    let (v, idx, round_state, counter) = vec_cipher_rj2::<F, 11>(c_rj2);

    // println!("{:?}", v);
    // println!("{:?}", idx);
    // println!("{:?}", round_state);

    let mut output: Vec<F> = Vec::new();

    for count in 0..counter{ 
        let round_indices: Vec<usize> = round_state.iter().enumerate().filter(|(_, &rs)| rs == count).map(|(index, _)| index).collect(); 

        let idxs: Vec<usize> = round_indices.iter().map(|&i| idx[i]).collect();

        let select: Vec<F> = round_indices.iter().map(|&i| v[i]).collect();
        let wits: Vec<F> = idxs.iter().map(|&i| wit_f[i]).collect();

        let ip = linalg::inner_product(&select, &wits); 
        if !haystack_r2j.contains(&ip) {
            println!("Count: {:?}", count);
            println!("{:?}", ip);
             println!("{:?}", idxs);
            println!("{:?}", select);
            println!("{:?}", wits);
        }
        output.push(linalg::inner_product(&select, &wits));
    }
    // println!("{:?}", output);


    assert!(output.into_iter().all(|x|(haystack_r2j.contains(&x))));
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
    aes_utils::rotate_right_inplace(&mut aux_m_col[0], 1);
    aes_utils::rotate_right_inplace(&mut aux_m_col[1], 2);
    aes_utils::rotate_right_inplace(&mut aux_m_col[2], 3);
    aes_utils::rotate_right_inplace(&mut aux_m_col[3], 3);

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

pub fn gcm_final_xor<F: Field, const R: usize>(dst: &mut [F], v: &[F], r: F, r2: F) {
    let reg = registry::aes_gcm_block_offsets::<R>();
    let mut v_pos = 0;

    //XOR over 16 byte messages
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
