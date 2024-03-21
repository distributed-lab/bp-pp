#![allow(non_snake_case)]

use std::ops::{Add, Mul};
use k256::{ProjectivePoint, Scalar};
use k256::elliptic_curve::rand_core::{CryptoRng, RngCore};
use merlin::Transcript;
use crate::range_proof::reciprocal::{Proof, ReciprocalRangeProof, Witness};

const G_VEC_CIRCUIT_SZ: usize = 16;
pub const G_VEC_FULL_SZ: usize = 16;
const H_VEC_CIRCUIT_SZ: usize = 26;
pub const H_VEC_FULL_SZ: usize = 32;


/// Represents public information for reciprocal range proof protocol for [0..2^64) range.
#[derive(Clone, Debug)]
pub struct U64RangeProof {
    /// Will be used for the value commitment as: `commitment = x*g + s*h_vec[0]`
    pub g: ProjectivePoint,
    /// Dimension: `16`
    pub g_vec: Vec<ProjectivePoint>,

    /// Will be used for the value commitment as: `commitment = x*g + s*h_vec[0]`
    /// Dimension: `26+6=32`
    pub h_vec: Vec<ProjectivePoint>,
}

impl U64RangeProof {
    /// Count of digits of u64 in hex representation.
    const DIM_ND: usize = 16;
    /// Base (hex)
    const DIM_NP: usize = 16;

    /// Creates commitment for the private value and blinding: `commitment = x*g + s*h_vec[0]`
    pub fn commit_value(&self, x: u64, s: &Scalar) -> ProjectivePoint {
        self.g.mul(&Scalar::from(x)).add(&self.h_vec[0].mul(s))
    }

    /// Verifies that committed value in `v` lies in range [0..2^64).
    pub fn verify(&self, v: &ProjectivePoint, proof: Proof, t: &mut Transcript) -> bool {
        let reciprocal = ReciprocalRangeProof {
            dim_nd: Self::DIM_ND,
            dim_np: Self::DIM_NP,
            g: self.g,
            g_vec: self.g_vec.clone(),
            h_vec: self.h_vec[..H_VEC_CIRCUIT_SZ].to_vec(),
            g_vec_: vec![],
            h_vec_: self.h_vec[H_VEC_CIRCUIT_SZ..].to_vec(),
        };

        reciprocal.verify(v, proof, t)
    }

    /// Creates proof that values `x` with blinding `s` lies in [0..2^64).
    pub fn prove<R>(&self, x: u64, s: &Scalar, t: &mut Transcript, rng: &mut R) -> Proof
        where
            R: RngCore + CryptoRng
    {
        let digits = Self::u64_to_hex(x);
        let poles = Self::u64_to_hex_mapped(x);

        let reciprocal = ReciprocalRangeProof {
            dim_nd: Self::DIM_ND,
            dim_np: Self::DIM_NP,
            g: self.g,
            g_vec: self.g_vec.clone(),
            h_vec: self.h_vec[..H_VEC_CIRCUIT_SZ].to_vec(),
            g_vec_: vec![],
            h_vec_: self.h_vec[H_VEC_CIRCUIT_SZ..].to_vec(),
        };

        let witness = Witness {
            x: Scalar::from(x),
            s: s.clone(),
            m: poles,
            digits,
        };

        reciprocal.prove(&reciprocal.commit_value(&witness.x, &witness.s), witness, t, rng)
    }

    pub fn u64_to_hex(mut x: u64) -> Vec<Scalar> {
        (0..16).map(|_| {
            let val = Scalar::from(x % 16);
            x = x / 16;
            val
        }).collect::<Vec<Scalar>>()
    }

    pub fn u64_to_hex_mapped(mut x: u64) -> Vec<Scalar> {
        let mut result = vec![Scalar::ZERO; 16];

        (0..16).for_each(|_| {
            let digit = (x % 16) as usize;
            result[digit] = result[digit].add(Scalar::ONE);
            x = x / 16;
        });

        result
    }
}