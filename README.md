# Bulletproofs++ implementation on Rust

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Abstract

Present Rust library contains the implementation of Bulletproofs++ weight norm linear argument protocol, arithmetic
circuit protocol and reciprocal range proof protocol. Also, contains the uint64 range proof protocol as a primary
use-case for reciprocal range proofs.

## Example of usage

```rust
pub fn main() {
    let mut rand = OsRng::default();

    let x = 123456u64; // private value to create proof for.
    let s = Scalar::generate_biased(&mut rand); // blinding value

    // Base points
    let g = k256::ProjectivePoint::random(&mut rand);
    let g_vec = (0..G_VEC_FULL_SZ).map(|_| k256::ProjectivePoint::random(&mut rand)).collect::<Vec<ProjectivePoint>>();
    let h_vec = (0..H_VEC_FULL_SZ).map(|_| k256::ProjectivePoint::random(&mut rand)).collect::<Vec<ProjectivePoint>>();

    let public = range_proof::u64_proof::U64RangeProof {
        g,
        g_vec,
        h_vec: h_vec[..H_VEC_CIRCUIT_SZ].to_vec(),
        h_vec_: h_vec[H_VEC_CIRCUIT_SZ..].to_vec(),
    };

    // transcript will be used for challenge generation - to move from interactive to non-interactive protocol.
    // transcript should be the same but new instance for prover and verifier. 
    let mut pt = Transcript::new(b"u64 range proof"); 
    let proof = public.prove(x, &s, &mut pt, &mut rand);

    // value commitment: `commitment = x*g + s*h_vec[0]`
    let commitment = public.commit_value(x, &s);

    let mut vt = Transcript::new(b"u64 range proof");
    assert!(public.verify(&commitment, proof, &mut vt));
}
```