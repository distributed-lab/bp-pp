use k256::{AffinePoint, ProjectivePoint, Scalar};

use super::Proof;

/// Represent serializable version of arithmetic circuit proof (uses AffinePoint instead of ProjectivePoint).
#[derive(serde::Serialize, serde::Deserialize, Clone, Debug)]
pub struct SerializableProof {
    pub c_l: AffinePoint,
    pub c_r: AffinePoint,
    pub c_o: AffinePoint,
    pub c_s: AffinePoint,
    pub r: Vec<AffinePoint>,
    pub x: Vec<AffinePoint>,
    pub l: Vec<Scalar>,
    pub n: Vec<Scalar>,
}

impl From<&SerializableProof> for Proof {
    fn from(value: &SerializableProof) -> Self {
        return Proof {
            c_l: ProjectivePoint::from(&value.c_l),
            c_r: ProjectivePoint::from(&value.c_r),
            c_o: ProjectivePoint::from(&value.c_o),
            c_s: ProjectivePoint::from(&value.c_s),
            r: value.r.iter().map(|r_val| ProjectivePoint::from(r_val)).collect::<Vec<ProjectivePoint>>(),
            x: value.x.iter().map(|x_val| ProjectivePoint::from(x_val)).collect::<Vec<ProjectivePoint>>(),
            l: value.l.clone(),
            n: value.n.clone(),
        };
    }
}

impl From<&Proof> for SerializableProof {
    fn from(value: &Proof) -> Self {
        return SerializableProof {
            c_l: value.c_l.to_affine(),
            c_r: value.c_r.to_affine(),
            c_o: value.c_o.to_affine(),
            c_s: value.c_s.to_affine(),
            r: value.r.iter().map(|r_val| r_val.to_affine()).collect::<Vec<AffinePoint>>(),
            x: value.x.iter().map(|x_val| x_val.to_affine()).collect::<Vec<AffinePoint>>(),
            l: value.l.clone(),
            n: value.n.clone(),
        };
    }
}
