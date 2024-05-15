use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};

use k256::elliptic_curve::Group;
use k256::elliptic_curve::rand_core::OsRng;
use k256::ProjectivePoint;
use bp_pp::range_proof;
use bp_pp::range_proof::u64_proof::{G_VEC_FULL_SZ, H_VEC_FULL_SZ};

fn bench_range_proof(ct: &mut Criterion) {
    // Preparing common data
    let mut rand = OsRng::default();

    let x = 123456u64;
    let s = k256::Scalar::generate_biased(&mut rand);


    // Preparing Base points
    let g = k256::ProjectivePoint::random(&mut rand);
    let g_vec = (0..G_VEC_FULL_SZ).map(|_| k256::ProjectivePoint::random(&mut rand)).collect::<Vec<ProjectivePoint>>();
    let h_vec = (0..H_VEC_FULL_SZ).map(|_| k256::ProjectivePoint::random(&mut rand)).collect::<Vec<ProjectivePoint>>();

    let public = range_proof::u64_proof::U64RangeProofProtocol {
        g,
        g_vec,
        h_vec,
    };


    // Benching proofs generation
    ct.bench_function("prove", |b| {
        b.iter_batched(
            || {
                merlin::Transcript::new(b"u64 range proof")
            },
            |mut pt| public.prove(x, &s, &mut pt, &mut rand),
            BatchSize::SmallInput,
        )
    });


    // Preparing proof to bench verification
    let proof = public.prove(x, &s, &mut merlin::Transcript::new(b"u64 range proof"), &mut rand);
    let commitment = public.commit_value(x, &s);

    ct.bench_function("verify", |b| {
        b.iter_batched(
            || {
                (merlin::Transcript::new(b"u64 range proof"), proof.clone())
            },
            |(mut vt, proof)| assert!(public.verify(&commitment, proof, &mut vt)),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(benches, bench_range_proof);
criterion_main!(benches);