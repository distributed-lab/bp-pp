use std::ops::{Add, Mul};
use k256::Scalar;

pub fn reduce<T: Copy>(v: &Vec<T>) -> (Vec<T>, Vec<T>) {
    let mut res0 = Vec::new();
    let mut res1 = Vec::new();

    for i in 0..v.len() {
        if i % 2 == 0 {
            res0.push(v[i]);
        } else {
            res1.push(v[i]);
        }
    }

    (res0, res1)
}


pub fn weight_vector_mul< T: Copy + Mul< Scalar, Output = T> + Add<Output = T> >(a: &Vec<T>, b: &Vec<Scalar>, weight: &Scalar) -> T {
    let mut result= a[0].mul(b[0].mul(weight));
    let mut exp = weight.mul(weight);

    for i in 1..a.len() {
        result = result.add(a[i].mul(b[i].mul(&exp)));
        exp = exp.mul(weight);
    }

    return result;
}

pub fn vector_mul< T: Copy + Mul<Scalar, Output = T> + Add<Output = T>>(a: &Vec<T>, b: &Vec<Scalar>) -> T {
    let mut result = a[0].mul(b[0]);

    for i in 1..a.len() {
        result = result.add(a[i].mul(b[i]));
    }

    return result;
}


pub fn vector_mul_on_scalar<'a, T: Copy + Mul<&'a Scalar, Output = T>>(a: &Vec<T>, s: &'a Scalar) -> Vec<T> {
    let mut result = vec![];

    for i in 0..a.len() {
        result.push(a[i].mul(s));
    }

    return result;
}


pub fn vector_add<T: Copy + Add<Output = T>>(a: &Vec<T>, b: &Vec<T>) -> Vec<T> {
    let mut result = vec![];

    for i in 0..a.len() {
        result.push(a[i].add(b[i]));
    }

    return result;
}
