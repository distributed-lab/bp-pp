use k256::{AffinePoint, ProjectivePoint};

use crate::circuit;

use super::reciprocal::Proof;


/// Represent serializable version of reciprocal proof (uses AffinePoint instead of ProjectivePoint
/// and serialized version of circuit proof).
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct SerializableProof {
    pub circuit_proof: circuit::SerializableProof,
    pub r: AffinePoint,
}

impl From<&SerializableProof> for Proof {
    fn from(value: &SerializableProof) -> Self {
        return Proof {
            circuit_proof: circuit::Proof::from(&value.circuit_proof),
            r: ProjectivePoint::from(value.r),
        };
    }
}

impl From<&Proof> for SerializableProof {
    fn from(value: &Proof) -> Self {
        return SerializableProof {
            circuit_proof: circuit::SerializableProof::from(&value.circuit_proof),
            r: value.r.to_affine(),
        };
    }
}
