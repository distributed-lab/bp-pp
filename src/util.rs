use std::ops::{Add, Mul, Sub};
use k256::elliptic_curve::Field;
use k256::Scalar;
use sha2::digest::generic_array::arr::Inc;

pub fn reduce<T: Copy>(v: &Vec<T>) -> (Vec<T>, Vec<T>) {
    let res0 = v.iter().
        enumerate().
        filter(|(i, _)| *i as i32 % 2 == 0).
        map(|(_, x)| *x).
        collect::<Vec<T>>();


    let res1 = v.iter().
        enumerate().
        filter(|(i, _)| *i as i32 % 2 == 1).
        map(|(_, x)| *x).
        collect::<Vec<T>>();

    (res0, res1)
}


pub fn weight_vector_mul<T: Copy + Mul<Scalar, Output=T> + Add<Output=T>>(a: &Vec<T>, b: &Vec<Scalar>, weight: &Scalar) -> T {
    let mut exp = weight.clone();
    let mut result = a[0].mul(b[0].mul(exp));

    (1..a.len()).
        map(|i| {
            exp = exp.mul(weight);
            return a[i].mul(b[i].mul(&exp));
        }).
        collect::<Vec<T>>().
        iter().
        for_each(|x| result = result.add(*x));

    return result;
}

pub fn vector_mul<T: Copy + Mul<Scalar, Output=T> + Add<Output=T>>(a: &Vec<T>, b: &Vec<Scalar>) -> T {
    let mut result = a[0].mul(b[0]);

    (1..a.len()).
        map(|i| a[i].mul(b[i])).
        collect::<Vec<T>>().
        iter().
        for_each(|x| result = result.add(*x));

    return result;
}


pub fn vector_mul_on_scalar<'a, T: Copy + Mul<&'a Scalar, Output=T>>(a: &Vec<T>, s: &'a Scalar) -> Vec<T> {
    a.iter().map(|x| x.mul(s)).collect::<Vec<T>>()
}


pub fn vector_add<T: Copy + Add<Output=T>>(a: &Vec<T>, b: &Vec<T>) -> Vec<T> {
    (0..a.len()).map(|i| a[i].add(b[i])).collect::<Vec<T>>()
}

pub fn vector_sub<'a, T: Copy + Sub<Output=T>>(a: &'a Vec<T>, b: &'a Vec<T>) -> Vec<T> {
    (0..a.len()).map(|i| a[i].sub(b[i])).collect::<Vec<T>>()
}

pub fn e(v: &Scalar, n: usize) -> Vec<Scalar> {
    let mut buf = Scalar::ONE;

    (0..n).map(|i| {
        let val = buf;
        buf = buf.mul(v);
        val
    }).collect::<Vec<Scalar>>()
}

pub fn pow(s: &Scalar, n: usize) -> Scalar {
    return s.pow(&[n as u64]);
}

pub fn vector_hadamard_mul<T: Copy + Mul<Scalar, Output=T>>(a: &Vec<T>, b: &Vec<Scalar>) -> Vec<T> {
    a.iter().enumerate().map(|(i, x)| x.mul(b[i])).collect::<Vec<T>>()
}

pub fn vector_tensor_mul<'a, T: Copy + Mul<&'a Scalar, Output=T>>(a: &'a Vec<T>, b: &'a Vec<Scalar>) -> Vec<T> {
    b.iter().map(|x| vector_mul_on_scalar(&a, x)).collect::<Vec<Vec<T>>>().concat()
}

pub fn diag_inv(x: &Scalar, n: usize) -> Vec<Vec<Scalar>> {
    let x_inv = x.invert().unwrap();
    let mut val = Scalar::ONE;

    (0..n).map(|i|
        (0..n).map(|j|
            if i == j {
                val = val.mul(x_inv);
                val
            } else {
                Scalar::ZERO
            }
        ).collect::<Vec<Scalar>>()
    ).collect::<Vec<Vec<Scalar>>>()
}

pub fn vector_mul_on_matrix<T: Copy + Mul<Scalar, Output=T> + Add<Output=T>>(a: &Vec<T>, m: &Vec<Vec<Scalar>>) -> Vec<T> {
    m.iter().map(|v| vector_mul(a, v)).collect::<Vec<T>>()
}

pub fn minus<T: Copy + Mul<Scalar, Output=T>>(v: &T) -> T {
    v.mul(Scalar::ZERO.sub(&Scalar::ONE))
}