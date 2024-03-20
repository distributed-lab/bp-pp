use std::ops::{Add, Mul, Sub};
use k256::{FieldBytes, ProjectivePoint, Scalar};
use k256::elliptic_curve::group::GroupEncoding;
use k256::elliptic_curve::PrimeField;
use k256::elliptic_curve::rand_core::{CryptoRng, RngCore};
use merlin::Transcript;
use crate::util::*;
use crate::wnla;
use crate::wnla::WeightNormLinearArgument;

#[derive(Clone, Copy)]
pub enum PartitionType {
    LO,
    LL,
    LR,
    NO,
}

pub trait Partition {
    fn map(typ: PartitionType, index: usize) -> Option<usize>;
}

#[derive(Debug)]
pub struct Proof {
    pub c_l: ProjectivePoint,
    pub c_r: ProjectivePoint,
    pub c_o: ProjectivePoint,
    pub c_s: ProjectivePoint,
    pub r: Vec<ProjectivePoint>,
    pub x: Vec<ProjectivePoint>,
    pub l: Vec<Scalar>,
    pub n: Vec<Scalar>,
}

pub struct Witness {
    // k*dim_nv
    pub v: Vec<Vec<Scalar>>,
    // k
    pub s_v: Vec<Scalar>,
    // dim_nm
    pub w_l: Vec<Scalar>,
    // dim_nm
    pub w_r: Vec<Scalar>,
    // dim_no
    pub w_o: Vec<Scalar>,
}

pub struct ArithmeticCircuit {
    pub dim_nm: usize,
    pub dim_no: usize,
    pub k: usize,

    // dim_nl = dim_nv * k
    pub dim_nl: usize,
    // Count of witness vectors v.
    pub dim_nv: usize,
    // dim_nw = dim_nm + dim_nm + n_o
    pub dim_nw: usize,

    pub g: ProjectivePoint,

    // dim_nm
    pub g_vec: Vec<ProjectivePoint>,
    // n_v+9
    pub h_vec: Vec<ProjectivePoint>,

    // dim_nm * dim_nw
    pub W_m: Vec<Vec<Scalar>>,
    // dim_nl * dim_nw
    pub W_l: Vec<Vec<Scalar>>,

    // dim_nm
    pub a_m: Vec<Scalar>,
    // dim_nl
    pub a_l: Vec<Scalar>,

    pub f_l: bool,
    pub f_m: bool,

    // Vectors of points that will be used in WNLA protocol
    // 2^n - dim_nm
    pub g_vec_: Vec<ProjectivePoint>,
    // 2^n - (dim_nv+9)
    pub h_vec_: Vec<ProjectivePoint>,
}

impl ArithmeticCircuit {
    pub fn commit(&self, v: &Vec<Scalar>, s: &Scalar) -> ProjectivePoint {
        self.
            g.mul(v[0]).
            add(self.h_vec[0].mul(s)).
            add(vector_mul(&self.h_vec[9..].to_vec(), &v[1..].to_vec()))
    }

    pub fn verify<T: RngCore + CryptoRng, P: Partition>(&self, v: &Vec<ProjectivePoint>, t: &mut Transcript, proof: Proof) -> bool {
        t.append_message(b"cl", proof.c_l.to_bytes().as_slice());
        t.append_message(b"cr", proof.c_r.to_bytes().as_slice());
        t.append_message(b"co", proof.c_o.to_bytes().as_slice());
        v.iter().for_each(|v_val| t.append_message(b"v", v_val.to_bytes().as_slice()));

        let rho = ArithmeticCircuit::get_challenge(b"rho", t);
        let lambda = ArithmeticCircuit::get_challenge(b"lambda", t);
        let beta = ArithmeticCircuit::get_challenge(b"beta", t);
        let delta = ArithmeticCircuit::get_challenge(b"delta", t);

        let mu = rho.mul(rho);

        let lambda_vec = self.collect_lambda(&lambda, &mu);
        let mu_vec = vector_mul_on_scalar(&e(&mu, self.dim_nm), &mu);

        let (c_nL, c_nR, c_nO, c_lL, c_lR, c_lO) = self.collect_c::<P>(&lambda_vec, &mu_vec, &mu);

        let mut v_ = ProjectivePoint::IDENTITY;
        (0..self.k).
            for_each(|i| v_ = v_.add(v[i].mul(self.linear_comb_coef(i, &lambda, &mu))));
        v_ = v_.mul(Scalar::from(2u32));

        t.append_message(b"cs", proof.c_s.to_bytes().as_slice());

        let tau = ArithmeticCircuit::get_challenge(b"tau", t);
        let tau_inv = tau.invert().unwrap();
        let tau2 = tau.mul(&tau);
        let tau3 = tau2.mul(&tau);

        println!("tau = {:?}", tau);

        let delta_inv = delta.invert().unwrap();

        let mut pn_tau = vector_mul_on_scalar(&c_nO, &tau3.mul(&delta_inv));
        pn_tau = vector_sub(&pn_tau, &vector_mul_on_scalar(&c_nL, &tau2));
        pn_tau = vector_add(&pn_tau, &vector_mul_on_scalar(&c_nR, &tau));

        let ps_tau = weight_vector_mul(&pn_tau, &pn_tau, &mu).
            add(&vector_mul(&lambda_vec, &self.a_l).mul(&tau3).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&mu_vec, &self.a_m).mul(&tau3).mul(&Scalar::from(2u32)));

        let pt = self.g.mul(ps_tau).add(vector_mul(&self.g_vec, &pn_tau));

        let cr_tau = vec![
            Scalar::ONE,
            tau_inv.mul(beta),
            tau.mul(beta),
            tau2.mul(beta),
            tau3.mul(beta),
            tau.mul(tau3).mul(beta),
            tau2.mul(tau3).mul(beta),
            tau3.mul(tau3).mul(beta),
            tau3.mul(tau3).mul(tau).mul(beta),
        ];

        let c_l0 = self.collect_cl0(&lambda, &mu);

        let mut cl_tau = vector_mul_on_scalar(&c_lO, &tau3.mul(&delta_inv));
        cl_tau = vector_sub(&cl_tau, &vector_mul_on_scalar(&c_lL, &tau2));
        cl_tau = vector_add(&cl_tau, &vector_mul_on_scalar(&c_lR, &tau));
        cl_tau = vector_mul_on_scalar(&cl_tau, &Scalar::from(2u32));
        cl_tau = vector_sub(&cl_tau, &c_l0);

        let mut c = Vec::from([&cr_tau[..], &cl_tau[..]].concat());

        println!("c = {:?}", c);

        let commitment = pt.
            add(&proof.c_s.mul(&tau_inv)).
            sub(&proof.c_o.mul(&delta)).
            add(&proof.c_l.mul(&tau)).
            sub(&proof.c_r.mul(&tau2)).
            add(&v_.mul(&tau3));

        println!("com = {:?}", commitment.to_affine());

        while c.len() < self.h_vec.len() + self.h_vec_.len() {
            c.push(Scalar::ZERO);
        }

        let wnla = WeightNormLinearArgument {
            g: ProjectivePoint::from(self.g),
            g_vec: Vec::from([&self.g_vec[..], &self.g_vec_[..]].concat()),
            h_vec: Vec::from([&self.h_vec[..], &self.h_vec_[..]].concat()),
            c,
            rho,
            mu,
        };


        println!("mu = {:?}", mu);

        return wnla.verify(&commitment, t, wnla::Proof {
            r: proof.r,
            x: proof.x,
            l: proof.l,
            n: proof.n,
        });
    }

    pub fn prove<T: RngCore + CryptoRng, P: Partition>(&self, v: &Vec<ProjectivePoint>, witness: Witness, t: &mut Transcript, rng: &mut T) -> Proof {
        let ro = vec![
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::ZERO,
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::ZERO,
        ];

        let rl = vec![
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::ZERO,
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::ZERO,
            Scalar::ZERO,
        ];

        let rr = vec![
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::ZERO,
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::generate_biased(rng),
            Scalar::ZERO,
            Scalar::ZERO,
            Scalar::ZERO,
        ];

        let nl = witness.w_l;
        let nr = witness.w_r;

        let no = (0..self.dim_nm).map(|j|
            if let Some(i) = P::map(PartitionType::NO, j) {
                Scalar::from(witness.w_o[i])
            } else {
                Scalar::ZERO
            }
        ).collect::<Vec<Scalar>>();

        let lo = (0..self.dim_nv).map(|j|
            if let Some(i) = P::map(PartitionType::LO, j) {
                Scalar::from(witness.w_o[i])
            } else {
                Scalar::ZERO
            }
        ).collect::<Vec<Scalar>>();

        let ll = (0..self.dim_nv).map(|j|
            if let Some(i) = P::map(PartitionType::LL, j) {
                Scalar::from(witness.w_o[i])
            } else {
                Scalar::ZERO
            }
        ).collect::<Vec<Scalar>>();

        let lr = (0..self.dim_nv).map(|j|
            if let Some(i) = P::map(PartitionType::LR, j) {
                Scalar::from(witness.w_o[i])
            } else {
                Scalar::ZERO
            }
        ).collect::<Vec<Scalar>>();

        let co =
            vector_mul(&self.h_vec, &Vec::from([&ro[..], &lo[..]].concat())).
                add(vector_mul(&self.g_vec, &no));

        let cl =
            vector_mul(&self.h_vec, &Vec::from([&rl[..], &ll[..]].concat())).
                add(vector_mul(&self.g_vec, &nl));

        let cr =
            vector_mul(&self.h_vec, &Vec::from([&rr[..], &lr[..]].concat())).
                add(vector_mul(&self.g_vec, &nr));

        t.append_message(b"cl", cl.to_bytes().as_slice());
        t.append_message(b"cr", cr.to_bytes().as_slice());
        t.append_message(b"co", co.to_bytes().as_slice());
        v.iter().for_each(|v_val| t.append_message(b"v", v_val.to_bytes().as_slice()));

        let rho = ArithmeticCircuit::get_challenge(b"rho", t);
        let lambda = ArithmeticCircuit::get_challenge(b"lambda", t);
        let beta = ArithmeticCircuit::get_challenge(b"beta", t);
        let delta = ArithmeticCircuit::get_challenge(b"delta", t);

        let mu = rho.mul(rho);

        let lambda_vec = self.collect_lambda(&lambda, &mu);
        let mu_vec = vector_mul_on_scalar(&e(&mu, self.dim_nm), &mu);


        // Check M matrix calculated ok

        // let (M_lnL, M_mnL, M_lnR, M_mnR) = self.collect_m_rl();
        // let (M_lnO, M_mnO, M_llL, M_mlL, M_llR, M_mlR, M_llO, M_mlO) = self.collect_m_o::<P>();

        // let mut W_lw = vector_add(&matrix_mul_on_vector(&lo, &M_llO), &matrix_mul_on_vector(&no, &M_lnO));
        // W_lw = vector_add(&W_lw, &vector_add(&matrix_mul_on_vector(&ll, &M_llL), &matrix_mul_on_vector(&nl, &M_lnL)));
        // W_lw = vector_add(&W_lw, &vector_add(&matrix_mul_on_vector(&lr, &M_llR), &matrix_mul_on_vector(&nr, &M_lnR)));
        //println!("Wl*w = {:?}", W_lw);

        // let mut W_mw = vector_add(&matrix_mul_on_vector(&lo, &M_mlO), &matrix_mul_on_vector(&no, &M_mnO));
        // W_mw = vector_add(&W_mw, &vector_add(&matrix_mul_on_vector(&ll, &M_mlL), &matrix_mul_on_vector(&nl, &M_mnL)));
        // W_mw = vector_add(&W_mw, &vector_add(&matrix_mul_on_vector(&lr, &M_mlR), &matrix_mul_on_vector(&nr, &M_mnR)));
        // println!("Wm*w = {:?}", W_mw);


        let (c_nL, c_nR, c_nO, c_lL, c_lR, c_lO) = self.collect_c::<P>(&lambda_vec, &mu_vec, &mu);

        let ls = (0..self.dim_nv).map(|_| Scalar::generate_biased(rng)).collect();
        let ns = (0..self.dim_nm).map(|_| Scalar::generate_biased(rng)).collect();

        let mut v_0 = Scalar::ZERO;
        (0..self.k).
            for_each(|i| v_0 = v_0.add(witness.v[i][0].mul(self.linear_comb_coef(i, &lambda, &mu))));
        v_0 = v_0.mul(Scalar::from(2u32));

        let mut rv = vec![Scalar::ZERO; 9];
        (0..self.k).
            for_each(|i| rv[0] = rv[0].add(witness.s_v[i].mul(self.linear_comb_coef(i, &lambda, &mu))));
        rv[0] = rv[0].mul(Scalar::from(2u32));

        let mut v_1 = vec![Scalar::ZERO; self.dim_nv - 1];
        (0..self.k).
            for_each(|i| v_1 = vector_add(&v_1, &vector_mul_on_scalar(&witness.v[i][1..].to_vec(), &self.linear_comb_coef(i, &lambda, &mu))));
        v_1 = vector_mul_on_scalar(&v_1, &Scalar::from(2u32));

        let c_l0 = self.collect_cl0(&lambda, &mu);

        // [-2 -1 0 1 2 4 5 6] -> f(tau) coefficients vector
        let mut f_ = vec![Scalar::ZERO; 8];

        let delta2 = delta.mul(&delta);
        let delta_inv = delta.invert().unwrap();

        // -2
        f_[0] = minus(&weight_vector_mul(&ns, &ns, &mu));

        // -1
        f_[1] = vector_mul(&c_l0, &ls).
            add(delta.mul(&Scalar::from(2u32)).mul(&weight_vector_mul(&ns, &no, &mu)));

        // 0
        f_[2] = minus(&vector_mul(&c_lR, &ls).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&c_l0, &lo).mul(&delta)).
            sub(&weight_vector_mul(&ns, &vector_add(&nl, &c_nR), &mu).mul(Scalar::from(2u32))).
            sub(&weight_vector_mul(&no, &no, &mu).mul(&delta2));

        //1
        f_[3] = vector_mul(&c_lL, &ls).mul(&Scalar::from(2u32)).
            add(&vector_mul(&c_lR, &lo).mul(&delta).mul(&Scalar::from(2u32))).
            add(&vector_mul(&c_l0, &ll)).
            add(&weight_vector_mul(&ns, &vector_add(&nr, &c_nL), &mu).mul(&Scalar::from(2u32))).
            add(&weight_vector_mul(&no, &vector_add(&nl, &c_nR), &mu).mul(&Scalar::from(2u32)).mul(&delta));

        // 2
        f_[4] = weight_vector_mul(&c_nR, &c_nR, &mu).
            sub(&vector_mul(&c_lO, &ls).mul(&delta_inv).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&c_lL, &lo).mul(&delta).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&c_lR, &ll).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&c_l0, &lr)).
            sub(&weight_vector_mul(&ns, &c_nO, &mu).mul(&delta_inv).mul(&Scalar::from(2u32))).
            sub(&weight_vector_mul(&no, &vector_add(&nr, &c_nL), &mu).mul(&delta).mul(&Scalar::from(2u32))).
            sub(&weight_vector_mul(&vector_add(&nl, &c_nR), &vector_add(&nl, &c_nR), &mu));

        // 3 should be zero

        // 4
        f_[5] = weight_vector_mul(&c_nO, &c_nR, &mu).mul(&delta_inv).mul(&Scalar::from(2u32)).
            add(&weight_vector_mul(&c_nL, &c_nL, &mu)).
            sub(&vector_mul(&c_lO, &ll).mul(&delta_inv).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&c_lL, &lr).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&c_lR, &v_1).mul(&Scalar::from(2u32))).
            sub(&weight_vector_mul(&vector_add(&nl, &c_nR), &c_nO, &mu).mul(&delta_inv).mul(&Scalar::from(2u32))).
            sub(&weight_vector_mul(&vector_add(&nr, &c_nL), &vector_add(&nr, &c_nL), &mu));

        // 5
        f_[6] = minus(&weight_vector_mul(&c_nO, &c_nL, &mu).mul(&delta_inv).mul(&Scalar::from(2u32))).
            add(&vector_mul(&c_nO, &lr).mul(&delta_inv).mul(&Scalar::from(2u32))).
            add(&vector_mul(&c_lL, &v_1).mul(&Scalar::from(2u32))).
            add(&weight_vector_mul(&vector_add(&nr, &c_nL), &c_nO, &mu).mul(&delta_inv).mul(&Scalar::from(2u32)));

        // 6
        f_[7] = minus(&vector_mul(&c_lO, &v_1).mul(&delta_inv).mul(&Scalar::from(2u32)));


        // TODO remove
        let f3 = vector_mul(&lambda_vec, &self.a_l).sub(&vector_mul(&mu_vec, &self.a_m)).mul(&Scalar::from(2u32)).
            sub(&weight_vector_mul(&c_nL, &c_nR, &mu).mul(&Scalar::from(2u32))).
            add(&v_0).
            add(&vector_mul(&c_lO, &lo).mul(&Scalar::from(2u32))).
            add(&vector_mul(&c_lL, &ll).mul(&Scalar::from(2u32))).
            add(&vector_mul(&c_lR, &lr).mul(&Scalar::from(2u32))).
            add(&vector_mul(&c_l0, &v_1)).
            add(&weight_vector_mul(&no, &c_nO, &mu).mul(&Scalar::from(2u32))).
            add(&weight_vector_mul(&vector_add(&nl, &c_nR), &vector_add(&nr, &c_nL), &mu).mul(&Scalar::from(2u32)));

        println!("f t^3 coef = {:?}", f3);

        let beta_inv = beta.invert().unwrap();

        let rs = vec![
            f_[1].add(ro[1].mul(&delta).mul(&beta)),
            f_[0].mul(&beta_inv),
            ro[0].mul(&delta).add(&f_[2]).mul(&beta_inv).sub(&rl[1]),
            f_[3].sub(&rl[0]).mul(&beta_inv).add(&ro[2].mul(&delta).add(rr[1])),
            f_[4].add(&rr[0]).mul(&beta_inv).add(&ro[3].mul(&delta).sub(rl[2])),
            minus(&rv[0].mul(&beta_inv)),
            f_[5].mul(&beta_inv).add(&ro[5].mul(&delta)).add(&rr[3]).sub(&rl[4]),
            f_[6].mul(&beta_inv).add(&rr[4]).add(&ro[6].mul(&delta)).sub(&rl[5]),
            f_[7].mul(&beta_inv).add(&ro[7].mul(&delta)).sub(&rl[6]).add(&rr[5]),
        ];

        let cs = vector_mul(&self.h_vec, &Vec::from([&rs[..], &ls[..]].concat())).
            add(vector_mul(&self.g_vec, &ns));

        t.append_message(b"cs", cs.to_bytes().as_slice());

        let tau = ArithmeticCircuit::get_challenge(b"tau", t);
        let tau_inv = tau.invert().unwrap();
        let tau2 = tau.mul(&tau);
        let tau3 = tau2.mul(&tau);

        println!("tau = {:?}", tau);

        let mut l = vector_mul_on_scalar(&Vec::from([&rs[..], &ls[..]].concat()), &tau_inv);
        l = vector_sub(&l, &vector_mul_on_scalar(&Vec::from([&ro[..], &lo[..]].concat()), &delta));
        l = vector_add(&l, &vector_mul_on_scalar(&Vec::from([&rl[..], &ll[..]].concat()), &tau));
        l = vector_sub(&l, &vector_mul_on_scalar(&Vec::from([&rr[..], &lr[..]].concat()), &tau2));
        l = vector_add(&l, &vector_mul_on_scalar(&Vec::from([&rv[..], &v_1[..]].concat()), &tau3));

        let mut pn_tau = vector_mul_on_scalar(&c_nO, &tau3.mul(&delta_inv));
        pn_tau = vector_sub(&pn_tau, &vector_mul_on_scalar(&c_nL, &tau2));
        pn_tau = vector_add(&pn_tau, &vector_mul_on_scalar(&c_nR, &tau));

        let ps_tau = weight_vector_mul(&pn_tau, &pn_tau, &mu).
            add(&vector_mul(&lambda_vec, &self.a_l).mul(&tau3).mul(&Scalar::from(2u32))).
            sub(&vector_mul(&mu_vec, &self.a_m).mul(&tau3).mul(&Scalar::from(2u32)));

        let mut n_tau = vector_mul_on_scalar(&ns, &tau_inv);
        n_tau = vector_sub(&n_tau, &vector_mul_on_scalar(&no, &delta));
        n_tau = vector_add(&n_tau, &vector_mul_on_scalar(&nl, &tau));
        n_tau = vector_sub(&n_tau, &vector_mul_on_scalar(&nr, &tau2));

        let mut n = vector_add(&pn_tau, &n_tau);

        let cr_tau = vec![
            Scalar::ONE,
            tau_inv.mul(beta),
            tau.mul(beta),
            tau2.mul(beta),
            tau3.mul(beta),
            tau.mul(tau3).mul(beta),
            tau2.mul(tau3).mul(beta),
            tau3.mul(tau3).mul(beta),
            tau3.mul(tau3).mul(tau).mul(beta),
        ];

        let mut cl_tau = vector_mul_on_scalar(&c_lO, &tau3.mul(&delta_inv));
        cl_tau = vector_sub(&cl_tau, &vector_mul_on_scalar(&c_lL, &tau2));
        cl_tau = vector_add(&cl_tau, &vector_mul_on_scalar(&c_lR, &tau));
        cl_tau = vector_mul_on_scalar(&cl_tau, &Scalar::from(2u32));
        cl_tau = vector_sub(&cl_tau, &c_l0);

        let mut c = Vec::from([&cr_tau[..], &cl_tau[..]].concat());
        println!("c = {:?}", c);

        let v = ps_tau.add(&tau3.mul(&v_0));

        let commitment = self.g.mul(v).
            add(&vector_mul(&self.h_vec, &l)).
            add(&vector_mul(&self.g_vec, &n));

        //c = [Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000001)), Scalar(Uint(0xC047368FAF68B8D7AEDC273588B38314971615A7DFAD631B90D091DCB35AC4A0)), Scalar(Uint(0x2A002C66B95660DA7CB97C452388ADA02D087F8CDBE9A5FC327FA5D1C0BE2885)), Scalar(Uint(0x1AFDDF56010588FADFBC9103FAF2F9CBE826A4F8E267F8A82CCA0840D022E147)), Scalar(Uint(0x78B0B00A678065E8A29FBB0935E9C80343F7AD4514DA341CAE59931009093E03)), Scalar(Uint(0xBA4DE638473052E6DFDC510347DFBF7D36A7B15D309B692900839DE59CD33DFA)), Scalar(Uint(0xB786E94A929AE9606DE7E886914C614708F974BF1B2BF33C5BAE6B6D4899087E)), Scalar(Uint(0x8C914AAB27774A20F4C36B00BED5245758A06409FB635FD71A1EF1D2ABD9EB07)), Scalar(Uint(0xD65574439D89237B8B98DD2C6FB806E6D6D1F8044C11F5A2466B241867953124)), Scalar(Uint(0x927C73DC6EC9D0DC8A1304A8601AE610DF02717969BF2C7914825512F38CFE61)), Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000000))]
        //c = [Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000001)), Scalar(Uint(0xC047368FAF68B8D7AEDC273588B38314971615A7DFAD631B90D091DCB35AC4A0)), Scalar(Uint(0x2A002C66B95660DA7CB97C452388ADA02D087F8CDBE9A5FC327FA5D1C0BE2885)), Scalar(Uint(0x1AFDDF56010588FADFBC9103FAF2F9CBE826A4F8E267F8A82CCA0840D022E147)), Scalar(Uint(0x78B0B00A678065E8A29FBB0935E9C80343F7AD4514DA341CAE59931009093E03)), Scalar(Uint(0xBA4DE638473052E6DFDC510347DFBF7D36A7B15D309B692900839DE59CD33DFA)), Scalar(Uint(0xB786E94A929AE9606DE7E886914C614708F974BF1B2BF33C5BAE6B6D4899087E)), Scalar(Uint(0x8C914AAB27774A20F4C36B00BED5245758A06409FB635FD71A1EF1D2ABD9EB07)), Scalar(Uint(0xD65574439D89237B8B98DD2C6FB806E6D6D1F8044C11F5A2466B241867953124)), Scalar(Uint(0x927C73DC6EC9D0DC8A1304A8601AE610DF02717969BF2C7914825512F38CFE61)), Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000000)), Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000000)), Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000000)), Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000000)), Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000000)), Scalar(Uint(0x0000000000000000000000000000000000000000000000000000000000000000))]

        println!("com = {:?}", commitment.to_affine());


        println!("Should be WNLA secret {:?}", v);

        while l.len() < self.h_vec.len() + self.h_vec_.len() {
            l.push(Scalar::ZERO);
            c.push(Scalar::ZERO);
        }

        while n.len() < self.g_vec.len() + self.g_vec_.len() {
            n.push(Scalar::ZERO);
        }

        let wnla = WeightNormLinearArgument {
            g: ProjectivePoint::from(self.g),
            g_vec: Vec::from([&self.g_vec[..], &self.g_vec_[..]].concat()),
            h_vec: Vec::from([&self.h_vec[..], &self.h_vec_[..]].concat()),
            c,
            rho,
            mu,
        };

        println!("com = {:?}", wnla.commit(&l, &n).to_affine());

        println!("mu = {:?}", mu);

        let proof_wnla = wnla.prove(&commitment, t, l, n);

        return Proof {
            c_l: cl,
            c_r: cr,
            c_o: co,
            c_s: cs,
            r: proof_wnla.r,
            x: proof_wnla.x,
            l: proof_wnla.l,
            n: proof_wnla.n,
        };
    }

    fn get_challenge(label: &'static [u8], t: &mut Transcript) -> Scalar {
        let mut buf = [0u8; 32];
        t.challenge_bytes(label, &mut buf);
        Scalar::from_repr(*FieldBytes::from_slice(&buf)).unwrap()
    }

    fn linear_comb_coef(&self, i: usize, lambda: &Scalar, mu: &Scalar) -> Scalar {
        let mut coef = Scalar::ZERO;
        if self.f_l {
            coef = coef.add(pow(&lambda, self.dim_nv * i))
        }

        if self.f_m {
            coef = coef.add(pow(&mu, self.dim_nv * i + 1))
        }

        coef
    }

    fn collect_cl0(&self, lambda: &Scalar, mu: &Scalar) -> Vec<Scalar> {
        let mut c_l0 = vec![Scalar::ZERO; self.dim_nv - 1];
        if self.f_l {
            c_l0 = e(&lambda, self.dim_nv)[1..].to_vec();
        }
        if self.f_m {
            c_l0 = vector_sub(&c_l0, &vector_mul_on_scalar(&e(&mu, self.dim_nv)[1..].to_vec(), &mu));
        }

        c_l0
    }

    fn collect_c<P: Partition>(&self, lambda_vec: &Vec<Scalar>, mu_vec: &Vec<Scalar>, mu: &Scalar) -> (Vec<Scalar>, Vec<Scalar>, Vec<Scalar>, Vec<Scalar>, Vec<Scalar>, Vec<Scalar>) {
        let (M_lnL, M_mnL, M_lnR, M_mnR) = self.collect_m_rl();
        let (M_lnO, M_mnO, M_llL, M_mlL, M_llR, M_mlR, M_llO, M_mlO) = self.collect_m_o::<P>();

        let mu_diag_inv = diag_inv(&mu, self.dim_nm);

        let c_nL = vector_mul_on_matrix(&vector_sub(&vector_mul_on_matrix(&lambda_vec, &M_lnL), &vector_mul_on_matrix(&mu_vec, &M_mnL)), &mu_diag_inv);
        let c_nR = vector_mul_on_matrix(&vector_sub(&vector_mul_on_matrix(&lambda_vec, &M_lnR), &vector_mul_on_matrix(&mu_vec, &M_mnR)), &mu_diag_inv);
        let c_nO = vector_mul_on_matrix(&vector_sub(&vector_mul_on_matrix(&lambda_vec, &M_lnO), &vector_mul_on_matrix(&mu_vec, &M_mnO)), &mu_diag_inv);

        let c_lL = vector_sub(&vector_mul_on_matrix(&lambda_vec, &M_llL), &vector_mul_on_matrix(&mu_vec, &M_mlL));
        let c_lR = vector_sub(&vector_mul_on_matrix(&lambda_vec, &M_llR), &vector_mul_on_matrix(&mu_vec, &M_mlR));
        let c_lO = vector_sub(&vector_mul_on_matrix(&lambda_vec, &M_llO), &vector_mul_on_matrix(&mu_vec, &M_mlO));

        return (c_nL, c_nR, c_nO, c_lL, c_lR, c_lO);
    }

    fn collect_lambda(&self, lambda: &Scalar, mu: &Scalar) -> Vec<Scalar> {
        let mut lambda_vec = e(lambda, self.dim_nl);
        if self.f_l && self.f_m {
            lambda_vec = vector_sub(
                &lambda_vec,
                &vector_add(
                    &vector_tensor_mul(&vector_mul_on_scalar(&e(lambda, self.dim_nv), mu), &e(&pow(mu, self.dim_nv), self.k)),
                    &vector_tensor_mul(&e(mu, self.dim_nv), &e(&pow(lambda, self.dim_nv), self.k)),
                ),
            );
        }

        return lambda_vec;
    }

    fn collect_m_rl(&self) -> (Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>) {
        let M_lnL = (0..self.dim_nl).map(|i| Vec::from(&self.W_l[i][..self.dim_nm])).collect::<Vec<Vec<Scalar>>>();
        let M_mnL = (0..self.dim_nm).map(|i| Vec::from(&self.W_m[i][..self.dim_nm])).collect::<Vec<Vec<Scalar>>>();
        let M_lnR = (0..self.dim_nl).map(|i| Vec::from(&self.W_l[i][self.dim_nm..self.dim_nm * 2])).collect::<Vec<Vec<Scalar>>>();
        let M_mnR = (0..self.dim_nm).map(|i| Vec::from(&self.W_m[i][self.dim_nm..self.dim_nm * 2])).collect::<Vec<Vec<Scalar>>>();
        return (M_lnL, M_mnL, M_lnR, M_mnR);
    }

    fn collect_m_o<P: Partition>(&self) -> (Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>, Vec<Vec<Scalar>>) {
        let W_lO = (0..self.dim_nl).map(|i| Vec::from(&self.W_l[i][self.dim_nm * 2..])).collect::<Vec<Vec<Scalar>>>();
        let W_mO = (0..self.dim_nm).map(|i| Vec::from(&self.W_m[i][self.dim_nm * 2..])).collect::<Vec<Vec<Scalar>>>();

        let map_f = |isz: usize, jsz: usize, typ: PartitionType, W_x: &Vec<Vec<Scalar>>| -> Vec<Vec<Scalar>>{
            (0..isz).map(|i|
                (0..jsz).map(|j|
                    if let Some(j_) = P::map(typ, j) {
                        Scalar::from(W_x[i][j_])
                    } else {
                        Scalar::ZERO
                    }
                ).collect::<Vec<Scalar>>()
            ).collect::<Vec<Vec<Scalar>>>()
        };

        let M_lnO = map_f(self.dim_nl, self.dim_nm, PartitionType::NO, &W_lO);
        let M_llL = map_f(self.dim_nl, self.dim_nv, PartitionType::LL, &W_lO);
        let M_llR = map_f(self.dim_nl, self.dim_nv, PartitionType::LR, &W_lO);
        let M_llO = map_f(self.dim_nl, self.dim_nv, PartitionType::LO, &W_lO);


        let M_mnO = map_f(self.dim_nm, self.dim_nm, PartitionType::NO, &W_mO);
        let M_mlL = map_f(self.dim_nm, self.dim_nv, PartitionType::LL, &W_mO);
        let M_mlR = map_f(self.dim_nm, self.dim_nv, PartitionType::LR, &W_mO);
        let M_mlO = map_f(self.dim_nm, self.dim_nv, PartitionType::LO, &W_mO);


        (M_lnO, M_mnO, M_llL, M_mlL, M_llR, M_mlR, M_llO, M_mlO)
    }
}
