use std::result;
use crate::aes_utils::{shiftrows, sbox, mixcolumns, xor, RC, rotate_right_inplace, RJ2};
use crate::aes_ks::{keyschedule,aes128_keyschedule, aes256_keyschedule};

use hex_literal::hex;

fn aes_round(mut state: [u8; 16], round_key: [u8; 16]) -> [u8; 16] {
    // Note: shiftrows before sbox
    state = shiftrows(state);
    state = sbox(state);
    state = mixcolumns(state);
    xor(state, round_key)
}

/// Naive implementation of AES-128
fn aes<const R: usize, const N: usize, const L: usize>(
    message: [u8; 16],
    key: [u8; L],
) -> [u8; 16] {
    debug_assert!((R == 11 && N == 4 && L == 16) || (R == 15 && N == 8 && L == 32));
    let keys = keyschedule::<R, N>(&key);

    let mut state = xor(message, keys[0]);
    for i in 1..R - 1 {
        state = aes_round(state, keys[i]);
    }
    let state = sbox(shiftrows(state));
    xor(state, keys[R - 1])
}

/// Naive implementation of AES-128
#[inline]
pub fn aes128(message: [u8; 16], key: [u8; 16]) -> [u8; 16] {
    aes::<11, 4, 16>(message, key)
}

/// Naive implementation of AES-256
#[inline]
pub fn aes256(message: [u8; 16], key: [u8; 32]) -> [u8; 16] {
    aes::<15, 8, 32>(message, key)
}

#[derive(Default)]
pub struct RoundTrace {
    // be careful here: SBOX is applied after shiftrows
    pub s_box: [u8; 16],
    pub m_col: [[u8; 16]; 5],
    pub start: [u8; 16],

    pub _s_row: [u8; 16],
    pub _aux_m_col: [[u8; 16]; 4],
}

/// The AES witness containing the full computation trace.
///
/// To have a rough idea of the sizes:
///
/// k_sch_s_box: 44
/// start: 160
/// final_s_box: 1s6
/// k_sch: 44 * 5
/// m_col: 144 * 5
#[derive(Default, Clone)]
pub struct AesCipherTrace {
    pub message: [u8; 16],
    pub key: [u8; 16],
    // cipher variables
    pub start: Vec<u8>,
    pub s_box: Vec<u8>,
    pub m_col: [Vec<u8>; 5],
    // last round
    pub output: [u8; 16],
    // key schedule permutations
    pub _keys: Vec<[u8; 16]>,
    // cipher permutations
    pub _s_row: Vec<u8>,
    pub _aux_m_col: [Vec<u8>; 4],
}

impl AesCipherTrace {
    pub fn add_round(&mut self, round_trace: &RoundTrace) {
        self._s_row.extend(&round_trace._s_row);
        self.s_box.extend(&round_trace.s_box);
        (0..5).for_each(|i| self.m_col[i].extend(&round_trace.m_col[i]));
        self.start.extend(&round_trace.start);
        (0..4).for_each(|i| self._aux_m_col[i].extend(&round_trace._aux_m_col[i]));
    }

    pub fn add_finalround(&mut self, trace: [[u8; 16]; 3]) {
        let [_final_s_row, final_s_box, output] = trace;
        self._s_row.extend(_final_s_row);
        self.s_box.extend(final_s_box);
        self.output = output;
    }

    pub fn new_aes128(message: [u8; 16], key: [u8; 16]) -> AesCipherTrace {
        let round_keys = aes128_keyschedule(&key);
        aes_trace(message, &round_keys)
    }

    pub fn new_aes256(message: [u8; 16], key: [u8; 32]) -> AesCipherTrace {
        let round_keys = aes256_keyschedule(&key);
        aes_trace(message, &round_keys)
    }
}

pub(crate) fn aes_trace<const R: usize>(
    message: [u8; 16],
    round_keys: &[[u8; 16]; R],
) -> AesCipherTrace {
    let mut witness = AesCipherTrace::default();
    witness.message = message;
    witness.key = round_keys[0];
    witness._keys = round_keys.to_vec();

    // first round: add key to message
    let mut round_state = xor(message, round_keys[0]);
    // put the first message into the key schedule artificially.
    // The verifier will do the same using the statement
    witness.start.extend(&round_state);
    #[allow(clippy::needless_range_loop)]
    for i in 1..R - 1 {
        let round_trace = aes_round_trace(round_state, round_keys[i]);
        witness.add_round(&round_trace);
        round_state = round_trace.start;
    }
    witness.add_finalround(final_round_trace(round_state, round_keys[R - 1]));
    witness
}

pub fn final_round_trace(state: [u8; 16], key: [u8; 16]) -> [[u8; 16]; 3] {
    let _s_row = shiftrows(state);
    let s_box = sbox(_s_row);
    let start = xor(s_box, key);
    [_s_row, s_box, start]
}

pub fn aes_round_trace(state: [u8; 16], key: [u8; 16]) -> RoundTrace {
    let mut trace = RoundTrace::default();
    // shiftrows
    trace._s_row = shiftrows(state);
    // sbox
    trace.s_box = sbox(trace._s_row);
    for (i, &x) in trace.s_box.iter().enumerate() {
        trace.m_col[0][i] = RJ2[x as usize];
    }
    // mixcolumns: generate the rotations of the vectors to xor.
    trace._aux_m_col[0] = trace.s_box;
    rotate_right_inplace(&mut trace._aux_m_col[0], 1);
    trace._aux_m_col[1] = trace.s_box;
    rotate_right_inplace(&mut trace._aux_m_col[1], 2);
    trace._aux_m_col[2] = trace.s_box;
    rotate_right_inplace(&mut trace._aux_m_col[2], 3);
    trace._aux_m_col[3] = trace.m_col[0];
    rotate_right_inplace(&mut trace._aux_m_col[3], 3);
    // mixcolumns
    trace.m_col[1] = xor(trace.m_col[0], trace._aux_m_col[0]);
    trace.m_col[2] = xor(trace.m_col[1], trace._aux_m_col[1]);
    trace.m_col[3] = xor(trace.m_col[2], trace._aux_m_col[2]);
    trace.m_col[4] = xor(trace.m_col[3], trace._aux_m_col[3]);
    trace.start = xor(trace.m_col[4], key);
    trace
}


#[test]
fn test_aes_round_trace() {
    let state = *b"\xc8\x16w\xbc\x9bz\xc9;%\x02y\x92\xb0&\x19\x96";
    let expected = *b"\xc6/\xe1\t\xf7^\xed\xc3\xccy9]\x84\xf9\xcf]";
    let round_key = *b"^9\x0f}\xf7\xa6\x92\x96\xa7U=\xc1\n\xa3\x1fk";

    let got = aes_round_trace(state, round_key);
    let naive = aes_round(state, round_key);

    assert_eq!(naive, expected);
    assert_eq!(got.start, expected);
}

#[test]
fn test_aes128() {
    let message: [u8; 16] = *b"\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff";
    let key: [u8; 16] = *b"\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f";

    let keys = aes128_keyschedule(&key);
    let expected = *b"\x00\x10 0@P`p\x80\x90\xa0\xb0\xc0\xd0\xe0\xf0";
    let state = xor(message, keys[0]);
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let state = sbox(state);
    let expected = *b"c\xca\xb7\x04\tS\xd0Q\xcd`\xe0\xe7\xbap\xe1\x8c";
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let state = shiftrows(state);
    let expected = *b"cS\xe0\x8c\t`\xe1\x04\xcdp\xb7Q\xba\xca\xd0\xe7";
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let state = mixcolumns(state);
    let expected = *b"_rd\x15W\xf5\xbc\x92\xf7\xbe;)\x1d\xb9\xf9\x1a";
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let state = xor(state, keys[1]);
    let expected = *b"\x89\xd8\x10\xe8\x85Z\xceh-\x18C\xd8\xcb\x12\x8f\xe4";
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let state = sbox(state);
    let state = shiftrows(state);
    let expected = *b"\xa7\xbe\x1ai\x97\xads\x9b\xd8\xc9\xcaE\x1fa\x8ba";
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let state = mixcolumns(state);
    let expected = *b"\xff\x87\x96\x841\xd8jQdQQ\xfaw:\xd0\t";
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let state = xor(state, keys[2]);
    let expected = *b"I\x15Y\x8fU\xe5\xd7\xa0\xda\xca\x94\xfa\x1f\nc\xf7";
    assert_eq!(
        expected,
        state,
        "\n{}\n{}",
        hex::encode(expected),
        hex::encode(state)
    );

    let got = aes128(message, key);
    let expected: [u8; 16] = *b"i\xc4\xe0\xd8j{\x040\xd8\xcd\xb7\x80p\xb4\xc5Z";
    assert_eq!(got, expected);
    let witness = AesCipherTrace::new_aes128(message, key);
    assert_eq!(witness.output, expected);
}

#[test]
fn test_aes256() {
    let message = *b"\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff";
    let key = *b"\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\
               \x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f";
    let expected = *b"\x8e\xa2\xb7\xca\x51\x67\x45\xbf\xea\xfc\x49\x90\x4b\x49\x60\x89";

    let round_keys = aes256_keyschedule(&key);
    let state = xor(message, round_keys[0]);
    assert_eq!(
        state,
        *b"\x00\x10\x20\x30\x40\x50\x60\x70\x80\x90\xa0\xb0\xc0\xd0\xe0\xf0"
    );
    let state = aes_round(state, round_keys[1]);
    assert_eq!(
        state,
        *b"\x4f\x63\x76\x06\x43\xe0\xaa\x85\xef\xa7\x21\x32\x01\xa4\xe7\x05"
    );
    let state = shiftrows(sbox(state));
    assert_eq!(
        state,
        *b"\x84\xe1\xfd\x6b\x1a\x5c\x94\x6f\xdf\x49\x38\x97\x7c\xfb\xac\x23"
    );
    let state = mixcolumns(state);
    assert_eq!(
        state,
        *b"\xbd\x2a\x39\x5d\x2b\x6a\xc4\x38\xd1\x92\x44\x3e\x61\x5d\xa1\x95"
    );

    let got = aes256(message, key);
    assert_eq!(got, expected);
}

#[test]
fn test_aes128_wiring() {
    let message: [u8; 16] = *b"\x00\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff";
    let key: [u8; 16] = *b"\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f";

    let got = aes128(message, key);
    let expected: [u8; 16] = *b"i\xc4\xe0\xd8j{\x040\xd8\xcd\xb7\x80p\xb4\xc5Z";
    assert_eq!(got, expected);
    let witness = AesCipherTrace::new_aes128(message, key);
    assert_eq!(witness.output, expected);

    // check a bit more in-depth the trace.
    let round_keys = aes128_keyschedule(&key);
    let start0 = xor(message, round_keys[0]);
    assert_eq!(start0, witness.start[..16]);
}
