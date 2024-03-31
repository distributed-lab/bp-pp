use k256::{AffinePoint, ProjectivePoint, Scalar};

use super::Proof;

/// Represent serializable version of  weight norm linear argument proof (uses AffinePoint instead of ProjectivePoint).
#[derive(Clone, Debug)]
pub struct SerializableProof {
    pub r: Vec<AffinePoint>,
    pub x: Vec<AffinePoint>,
    pub l: Vec<Scalar>,
    pub n: Vec<Scalar>,
}

impl From<&SerializableProof> for Proof {
    fn from(value: &SerializableProof) -> Self {
        return Proof {
            r: value.r.iter().map(ProjectivePoint::from).collect::<Vec<ProjectivePoint>>(),
            x: value.x.iter().map(ProjectivePoint::from).collect::<Vec<ProjectivePoint>>(),
            l: value.l.clone(),
            n: value.n.clone(),
        };
    }
}

impl From<&Proof> for SerializableProof {
    fn from(value: &Proof) -> Self {
        return SerializableProof {
            r: value.r.iter().map(|r_val| r_val.to_affine()).collect::<Vec<AffinePoint>>(),
            x: value.x.iter().map(|x_val| x_val.to_affine()).collect::<Vec<AffinePoint>>(),
            l: value.l.clone(),
            n: value.n.clone(),
        };
    }
}
