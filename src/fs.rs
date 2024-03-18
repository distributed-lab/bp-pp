use k256::{AffinePoint, ProjectivePoint, Scalar};
use k256::elliptic_curve::group::GroupEncoding;
use k256::elliptic_curve::PrimeField;
use sha2::{Sha256, Digest};

pub trait FiatShamir {
    fn add_point(&mut self, p: &k256::ProjectivePoint) -> &mut Self;
    fn add_scalar(&mut self, s: &k256::Scalar) -> &mut Self;
    fn get_challenge(&mut self) -> k256::Scalar;
}

pub struct SHA2FS {
    pub state: Vec<u8>,
}

impl SHA2FS {
    pub fn new() -> Self {
        return SHA2FS {
            state: Vec::new(),
        };
    }
}

impl FiatShamir for SHA2FS {
    fn add_point(&mut self, p: &ProjectivePoint) -> &mut Self {
        self.state.append(&mut Vec::from(p.to_bytes().as_slice()));
        self
    }

    fn add_scalar(&mut self, s: &Scalar) -> &mut Self {
        self.state.append(&mut Vec::from(s.to_bytes().as_slice()));
        self
    }
    fn get_challenge(&mut self) -> Scalar {
        self.state = Vec::from(Sha256::digest(self.state.as_slice()).as_slice());
        let bytes = k256::FieldBytes::from_slice(self.state.as_slice());
        k256::Scalar::from_repr(*bytes).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fs_sha2_works() {
        let mut fs = SHA2FS::new();
        fs.add_scalar(&k256::Scalar::from(1 as u32));
        fs.add_scalar(&k256::Scalar::from(2 as u32));

        let challenge = fs.get_challenge();
        println!("{}", k256::U256::from(challenge));

        fs.add_scalar(&k256::Scalar::from(3 as u32));

        let challenge = fs.get_challenge();
        println!("{}", k256::U256::from(challenge));

        let challenge = fs.get_challenge();
        println!("{}", k256::U256::from(challenge));
    }
}
