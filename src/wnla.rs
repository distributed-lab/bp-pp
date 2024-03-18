use std::ops::{Add, Mul};
use k256::{FieldBytes, ProjectivePoint, Scalar};
use k256::elliptic_curve::group::GroupEncoding;
use k256::elliptic_curve::group::prime::PrimeCurveAffine;
use k256::elliptic_curve::PrimeField;
use merlin::Transcript;
use crate::util::*;

pub struct WNLAPublic {
    pub g: ProjectivePoint,
    pub g_vec: Vec<ProjectivePoint>,
    pub h_vec: Vec<ProjectivePoint>,
    pub c: Vec<Scalar>,
    pub rho: Scalar,
    pub mu: Scalar,
}

#[derive(Debug)]
pub struct WNLAProof {
    pub r: Vec<ProjectivePoint>,
    pub x: Vec<ProjectivePoint>,
    pub l: Vec<Scalar>,
    pub n: Vec<Scalar>,
}

impl WNLAPublic {
    pub fn commit(&self, l: &Vec<Scalar>, n: &Vec<Scalar>) -> ProjectivePoint {
        let v = vector_mul(&self.c, l).add(weight_vector_mul(n, n, &self.mu));
        self.
            g.mul(v).
            add(vector_mul(&self.h_vec, l)).
            add(vector_mul(&self.g_vec, n))
    }
}

pub fn verify(public: &WNLAPublic, com: &ProjectivePoint, t: &mut Transcript, proof: WNLAProof) -> bool {
    if proof.x.len() != proof.r.len() {
        return false;
    }

    if proof.x.len() == 0 {
        return com.eq(&public.commit(&proof.l, &proof.n));
    }

    let (c0, c1) = reduce(&public.c);
    let (g0, g1) = reduce(&public.g_vec);
    let (h0, h1) = reduce(&public.h_vec);

    t.append_message(b"com", com.to_bytes().as_slice());
    t.append_message(b"x", proof.x.last().unwrap().to_bytes().as_slice());
    t.append_message(b"r", proof.r.last().unwrap().to_bytes().as_slice());
    t.append_u64(b"l.sz", public.h_vec.len() as u64);
    t.append_u64(b"n.sz", public.g_vec.len() as u64);

    let mut buf = [0u8; 32];
    t.challenge_bytes(b"y", &mut buf);
    let y = k256::Scalar::from_repr(*FieldBytes::from_slice(&buf)).unwrap();

    let h_ = vector_add(&h0, &vector_mul_on_scalar(&h1, &y));
    let g_ = vector_add(&vector_mul_on_scalar(&g0, &public.rho), &vector_mul_on_scalar(&g1, &y));
    let c_ = vector_add(&c0, &vector_mul_on_scalar(&c1, &y));


    let com_ = com.
        add(&proof.x[0].mul(y)).
        add(&proof.r[0].mul(y.mul(y).sub(&Scalar::ONE)));

    let public_ = WNLAPublic {
        g: public.g,
        g_vec: g_,
        h_vec: h_,
        c: c_,
        rho: public.mu,
        mu: public.mu.mul(&public.mu),
    };

    let proof_ = WNLAProof {
        r: proof.r[..proof.r.len() - 1].to_vec(),
        x: proof.r[..proof.x.len() - 1].to_vec(),
        l: proof.l,
        n: proof.n,
    };

    return verify(&public_, &com_, t, proof_);
}

pub fn prove(public: &WNLAPublic, com: &ProjectivePoint, t: &mut Transcript, l: Vec<Scalar>, n: Vec<Scalar>) -> WNLAProof {
    if l.len() + n.len() < 6 {
        return WNLAProof {
            r: vec![],
            x: vec![],
            l,
            n,
        };
    }

    let rho_inv = public.rho.invert().unwrap();

    let (c0, c1) = reduce(&public.c);
    let (l0, l1) = reduce(&l);
    let (n0, n1) = reduce(&n);
    let (g0, g1) = reduce(&public.g_vec);
    let (h0, h1) = reduce(&public.h_vec);

    let mu2 = public.mu.mul(&public.mu);

    let vx = weight_vector_mul(&n0, &n1, &mu2).
        mul(&rho_inv.mul(&Scalar::from(2u32))).
        add(&vector_mul(&c0, &l1)).
        add(&vector_mul(&c1, &l0));

    let vr = weight_vector_mul(&n1, &n1, &mu2).add(&vector_mul(&c1, &l1));

    let x = public.g.mul(vx).
        add(&vector_mul(&h0, &l1)).
        add(&vector_mul(&h1, &l0)).
        add(&vector_mul(&g0, &vector_mul_on_scalar(&n1, &public.rho))).
        add(&vector_mul(&g1, &vector_mul_on_scalar(&n0, &rho_inv)));

    let r = public.g.mul(vr).
        add(vector_mul(&h1, &l1)).
        add(vector_mul(&g1, &n1));


    t.append_message(b"com", com.to_bytes().as_slice());
    t.append_message(b"x", x.to_bytes().as_slice());
    t.append_message(b"r", r.to_bytes().as_slice());
    t.append_u64(b"l.sz", l.len() as u64);
    t.append_u64(b"n.sz", n.len() as u64);

    let mut buf = [0u8; 32];
    t.challenge_bytes(b"y", &mut buf);
    let y = k256::Scalar::from_repr(*FieldBytes::from_slice(&buf)).unwrap();

    let h_ = vector_add(&h0, &vector_mul_on_scalar(&h1, &y));
    let g_ = vector_add(&vector_mul_on_scalar(&g0, &public.rho), &vector_mul_on_scalar(&g1, &y));
    let c_ = vector_add(&c0, &vector_mul_on_scalar(&c1, &y));

    let l_ = vector_add(&l0, &vector_mul_on_scalar(&l1, &y));
    let n_ = vector_add(&vector_mul_on_scalar(&n0, &rho_inv), &vector_mul_on_scalar(&n1, &y));

    let public_ = WNLAPublic {
        g: public.g,
        g_vec: g_,
        h_vec: h_,
        c: c_,
        rho: public.mu,
        mu: mu2,
    };

    let mut proof = prove(&public_, &public_.commit(&l_, &n_), t, l_, n_);
    proof.r.push(r);
    proof.x.push(x);
    proof
}


#[cfg(test)]
mod tests {
    use k256::elliptic_curve::Group;
    use k256::elliptic_curve::rand_core::OsRng;
    use k256::Scalar;
    use super::*;

    #[test]
    fn wnla_works() {
        const N: i32 = 4;

        let mut rand = OsRng::default();

        let g = k256::ProjectivePoint::random(&mut rand);

        let mut g_vec = vec![];
        for i in 0..N {
            g_vec.push(k256::ProjectivePoint::random(&mut rand));
        }

        let mut h_vec = vec![];
        for i in 0..N {
            h_vec.push(k256::ProjectivePoint::random(&mut rand));
        }

        let mut c = vec![];
        for i in 0..N {
            c.push(k256::Scalar::generate_biased(&mut rand));
        }

        let rho = k256::Scalar::generate_biased(&mut rand);

        let public = WNLAPublic {
            g,
            g_vec,
            h_vec,
            c,
            rho,
            mu: rho.mul(&rho),
        };

        let l = vec![Scalar::from(1 as u32), Scalar::from(2 as u32), Scalar::from(3 as u32), Scalar::from(4 as u32)];
        let n = vec![Scalar::from(8 as u32), Scalar::from(7 as u32), Scalar::from(6 as u32), Scalar::from(5 as u32)];

        let commit = public.commit(&l, &n);
        let mut pt = Transcript::new(b"wnla test");

        let proof = prove(&public, &commit, &mut pt, l, n);
        println!("{:?}", &proof);

        let mut vt = Transcript::new(b"wnla test");
        println!("{}", verify(&public, &commit, &mut vt, proof))
    }
}
