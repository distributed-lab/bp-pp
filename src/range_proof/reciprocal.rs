#![allow(non_snake_case)]

use std::ops::{Add, Mul};
use k256::{ProjectivePoint, Scalar};
use k256::elliptic_curve::rand_core::{CryptoRng, RngCore};
use merlin::Transcript;
use crate::util::*;
use crate::{circuit, transcript};
use crate::circuit::{ArithmeticCircuit, PartitionType};

/// Represents reciprocal range-proof protocol witness.
#[derive(Clone, Debug)]
pub struct Witness {
    /// Private value in range [0..dim_np^dim_nd).
    pub x: Scalar,
    /// Blinding value
    pub s: Scalar,
    /// Witness vector: value multiplicities: i-th element corresponds to the 'i-digit' multiplicity
    pub m: Vec<Scalar>,
    /// Witness vector: value digits.
    pub digits: Vec<Scalar>,
}

/// Represents reciprocal range-proof: zk-proof that committed value lies in [0..dim_np^dim_nd) range.
#[derive(Clone, Debug)]
pub struct Proof {
    pub circuit_proof: circuit::Proof,
    pub r: ProjectivePoint,
}

/// Represents public reciprocal range proof information. Using this information and challenge both
/// prover and verifier can derive the arithmetic circuit.
#[derive(Clone, Debug)]
pub struct ReciprocalRangeProof {
    /// Count of private proles (size of committed value). Equals to: `dim_nm`. Also, `dim_nv = 1 + dim_nd`.
    pub dim_nd: usize,
    /// Count of public poles (number system base). Equals to: `dim_no`.
    pub dim_np: usize,

    /// Will be used for the value commitment: `commitment = x*g + s*h_vec[0]`
    pub g: ProjectivePoint,

    /// Dimension: `dim_nm`
    pub g_vec: Vec<ProjectivePoint>,
    /// Will be used for the value commitment: `commitment = x*g + s*h_vec[0]`
    /// Dimension: `dim_nv+9`
    pub h_vec: Vec<ProjectivePoint>,

    /// Additional points to be used in WNLA.
    /// Dimension: `2^n - dim_nm`
    pub g_vec_: Vec<ProjectivePoint>,
    /// Dimension: `2^n - (dim_nv+9)`
    pub h_vec_: Vec<ProjectivePoint>,
}

impl ReciprocalRangeProof {
    /// Creates commitment for the private value and blinding: `commitment = x*g + s*h_vec[0]`
    pub fn commit_value(&self, x: &Scalar, s: &Scalar) -> ProjectivePoint {
        self.g.mul(x).add(&self.h_vec[0].mul(s))
    }

    /// Creates commitment for the reciprocals and blinding: `commitment = s*h_vec[0] + <r, h_vec[9:]>`
    pub fn commit_poles(&self, r: &Vec<Scalar>, s: &Scalar) -> ProjectivePoint {
        self.h_vec[0].mul(s).add(&vector_mul(&self.h_vec[9..].to_vec(), &r))
    }

    /// Verifies zk-proof that committed value lies in [0..dim_np^dim_nd) range.
    pub fn verify(&self, commitment: &ProjectivePoint, proof: Proof, t: &mut Transcript) -> bool {
        transcript::app_point(b"reciprocal_commitment", commitment, t);
        let e = transcript::get_challenge(b"reciprocal_challenge", t);

        let partition = |typ: PartitionType, index: usize| -> Option<usize>{
            if typ == PartitionType::LL && index < self.dim_np {
                Some(index)
            } else {
                None
            }
        };

        let circuit = self.make_circuit(e, partition);

        let circuit_commitment = commitment.add(&proof.r);

        return circuit.verify(&vec![circuit_commitment], t, proof.circuit_proof);
    }

    /// Creates zk-proof that committed value lies in [0..dim_np^dim_nd) range.
    pub fn prove<R>(&self, commitment: &ProjectivePoint, witness: Witness, t: &mut Transcript, rng: &mut R) -> Proof
        where
            R: RngCore + CryptoRng
    {
        transcript::app_point(b"reciprocal_commitment", commitment, t);
        let e = transcript::get_challenge(b"reciprocal_challenge", t);

        let r = (0..self.dim_nd).map(|i|
            witness.digits[i].add(&e).invert().unwrap()
        ).collect::<Vec<Scalar>>();

        let r_blind = Scalar::generate_biased(rng);
        let r_com = self.commit_poles(&r, &r_blind);

        let mut v = vec![witness.x];
        r.iter().for_each(|r_val| v.push(*r_val));

        let w_l = witness.digits;
        let w_r = r;
        let w_o = witness.m;

        let partition = |typ: PartitionType, index: usize| -> Option<usize>{
            if typ == PartitionType::LL && index < self.dim_np {
                Some(index)
            } else {
                None
            }
        };

        let circuit = self.make_circuit(e, partition);

        let circuit_witness = circuit::Witness {
            v: vec![v],
            s_v: vec![witness.s.add(r_blind)],
            w_l,
            w_r,
            w_o,
        };

        let circuit_commitment = circuit.commit(&circuit_witness.v[0], &circuit_witness.s_v[0]);
        return Proof {
            circuit_proof: circuit.prove::<R>(&vec![circuit_commitment], circuit_witness, t, rng),
            r: r_com,
        };
    }

    fn make_circuit<P>(&self, e: Scalar, partition: P) -> ArithmeticCircuit<P>
        where
            P: Fn(PartitionType, usize) -> Option<usize>
    {
        let dim_nm = self.dim_nd;
        let dim_no = self.dim_np;

        let dim_nv = self.dim_nd + 1;
        let dim_nl = dim_nv;
        let dim_nw = self.dim_nd * 2 + self.dim_np;

        let a_m = vec![Scalar::ONE; dim_nm];

        let mut W_m = vec![vec![Scalar::ZERO; dim_nw]; dim_nm];
        (0..dim_nm).for_each(|i| W_m[i][i + dim_nm] = minus(&e));

        let a_l = vec![Scalar::ZERO; dim_nl];
        let base = Scalar::from(self.dim_np as u32);

        let mut W_l = vec![vec![Scalar::ZERO; dim_nw]; dim_nl];

        // fill for v-part in w vector
        (0..dim_nm).for_each(|i| W_l[0][i] = minus(&pow(&base, i)));

        // fill for r-part in w vector
        (0..dim_nm).for_each(|i|
            (0..dim_nm).for_each(|j| W_l[i + 1][j + dim_nm] = Scalar::ONE)
        );

        (0..dim_nm).for_each(|i| W_l[i + 1][i + dim_nm] = Scalar::ZERO);

        (0..dim_nm).for_each(|i|
            (0..dim_no).for_each(|j|
                W_l[i + 1][j + 2 * dim_nm] = minus(&(e.add(Scalar::from(j as u32)).invert().unwrap()))
            )
        );

        ArithmeticCircuit {
            dim_nm,
            dim_no,
            k: 1,
            dim_nl,
            dim_nv,
            dim_nw,
            g: self.g.clone(),
            g_vec: self.g_vec.clone(),
            h_vec: self.h_vec.clone(),
            W_m,
            W_l,
            a_m,
            a_l,
            f_l: true,
            f_m: false,
            g_vec_: self.g_vec_.clone(),
            h_vec_: self.h_vec_.clone(),
            partition,
        }
    }
}


