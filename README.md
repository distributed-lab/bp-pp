# Bulletproofs++ implementation on Rust

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Pull Requests welcome](https://img.shields.io/badge/PRs-welcome-ff69b4.svg?style=flat-square)](https://github.com/distributed-lab/bp-pp/issues)
<a href="https://github.com/distributed-lab/bp-pp">
<img src="https://img.shields.io/github/stars/distributed-lab/bp-pp?style=social"/>
</a>

## Abstract

Present Rust library contains the implementation of Bulletproofs++ weight norm linear argument protocol, arithmetic
circuit protocol and reciprocal range proof protocol. Also, contains the uint64 range proof protocol as a primary
use-case for reciprocal range proofs.

Implemented protocol has 2G points advantage over existing BP and BP+ protocols on proving of one 64-bit value and this
advantage will increase for more values per proof.

| Protocol | G  | F |
|----------|----|---|
| BP       | 16 | 5 |
| BP+      | 15 | 3 |
| Our BP++ | 13 | 3 |

This implementation uses [Merlin transcript](https://doc.dalek.rs/merlin/index.html) for challenges generation as was
recommended by Bulletproofs protocol authors.

All `Proof` data models has corresponding `SerializeProof` models where [serde](https://serde.rs/) `Serialize`
and `Deserialize` was implemented.

## Example of usage

Use [tests](./src/tests.rs) to run the provided example: 

```rust
pub fn main() {
    let mut rand = OsRng::default();

    let x = 123456u64; // private value to create proof for.
    let s = k256::Scalar::generate_biased(&mut rand); // blinding value

    // Base points
    let g = k256::ProjectivePoint::random(&mut rand);
    let g_vec = (0..G_VEC_FULL_SZ).map(|_| k256::ProjectivePoint::random(&mut rand)).collect::<Vec<ProjectivePoint>>();
    let h_vec = (0..H_VEC_FULL_SZ).map(|_| k256::ProjectivePoint::random(&mut rand)).collect::<Vec<ProjectivePoint>>();

    let public = range_proof::u64_proof::U64RangeProofProtocol {
        g,
        g_vec,
        h_vec,
    };

    // transcript will be used for challenge generation - to move from interactive to non-interactive protocol.
    // transcript should be the new instance but with same label for prover and verifier. 
    let mut pt = merlin::Transcript::new(b"u64 range proof");
    let proof = public.prove(x, &s, &mut pt, &mut rand);

    // value commitment: `commitment = x*g + s*h_vec[0]`
    let commitment = public.commit_value(x, &s);

    println!("{}", serde_json::to_string_pretty(&reciprocal::SerializableProof::from(&proof)).unwrap());

    let mut vt = merlin::Transcript::new(b"u64 range proof");
    assert!(public.verify(&commitment, proof, &mut vt));
}
```