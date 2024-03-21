#![allow(non_snake_case)]

use std::ops::{Add, Mul};
use k256::{ProjectivePoint, Scalar};
use k256::elliptic_curve::rand_core::{CryptoRng, RngCore};
use merlin::Transcript;
use crate::range_proof::reciprocal::{Proof, ReciprocalRangeProof, Witness};

pub struct U64RangeProof {
    // Base points

    // g and h_vec[0] will be used for the value commitment: commitment = x*g + s*h_vec[0]
    pub g: ProjectivePoint,

    // len = 16
    pub g_vec: Vec<ProjectivePoint>,
    // len = 26
    pub h_vec: Vec<ProjectivePoint>,

    // Additional points

    // len = 6
    pub h_vec_: Vec<ProjectivePoint>,
}

impl U64RangeProof {
    const DIM_ND: usize = 16;
    const DIM_NP: usize = 16;
    pub fn commit_value(&self, x: u64, s: &Scalar) -> ProjectivePoint {
        self.g.mul(&Scalar::from(x)).add(&self.h_vec[0].mul(s))
    }

    pub fn verify(&self, v: &ProjectivePoint, proof: Proof, t: &mut Transcript) -> bool {
        let reciprocal = ReciprocalRangeProof {
            dim_nd: Self::DIM_ND,
            dim_np: Self::DIM_NP,
            g: self.g,
            g_vec: self.g_vec.clone(),
            h_vec: self.h_vec.clone(),
            g_vec_: vec![],
            h_vec_: self.h_vec_.clone(),
        };

        reciprocal.verify(v, proof, t)
    }

    pub fn prove<T: RngCore + CryptoRng>(&self, x: u64, s: &Scalar, t: &mut Transcript, rng: &mut T) -> Proof {
        let digits = Self::u64_to_hex(x);
        let poles = Self::u64_to_hex_mapped(x);

        let reciprocal = ReciprocalRangeProof {
            dim_nd: Self::DIM_ND,
            dim_np: Self::DIM_NP,
            g: self.g,
            g_vec: self.g_vec.clone(),
            h_vec: self.h_vec.clone(),
            g_vec_: vec![],
            h_vec_: self.h_vec_.clone(),
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