use k256::{FieldBytes, ProjectivePoint, Scalar};
use k256::elliptic_curve::group::GroupEncoding;
use k256::elliptic_curve::PrimeField;
use merlin::Transcript;


pub fn app_point(label: &'static [u8], p: &ProjectivePoint, t: &mut Transcript) {
    t.append_message(label, p.to_bytes().as_slice());
}

pub fn get_challenge(label: &'static [u8], t: &mut Transcript) -> Scalar {
    let mut buf = [0u8; 32];
    t.challenge_bytes(label, &mut buf);
    Scalar::from_repr(*FieldBytes::from_slice(&buf)).unwrap()
}
