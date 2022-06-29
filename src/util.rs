/**
Implement traits to emulate vector-like element-wise addition for Vec<f64>
*/
use std::ops::{Add, AddAssign};
use std::vec;

pub trait ElementWiseAdd<T, Rhs = Self> {
    type Output;
    fn add(&self, rhs: &Rhs) -> Self::Output;
}

pub trait ElementWiseAddAssign<T, Rhs = Self> {
    fn add_assign(&mut self, rhs: &Rhs);
}

impl<T> ElementWiseAdd<T> for Vec<T>
where
    T: Copy + Add,
    Vec<T>: FromIterator<<T as Add>::Output>,
{
    type Output = Self;
    fn add(&self, other: &Self) -> Self::Output {
        self.iter()
            .zip(other.iter())
            .map(|(a, b)| *a + *b)
            .collect()
    }
}

impl<T> ElementWiseAdd<T> for &[T]
where
    T: Copy + Add,
    Vec<T>: FromIterator<<T as Add>::Output>,
{
    type Output = Vec<T>;
    fn add(&self, other: &Self) -> Vec<T> {
        self.iter()
            .zip(other.iter())
            .map(|(a, b)| *a + *b)
            .collect()
    }
}

impl<T> ElementWiseAddAssign<T> for Vec<T>
where
    T: Copy + AddAssign,
{
    fn add_assign(&mut self, other: &Self) {
        self.iter_mut()
            .zip(other.iter())
            .for_each(|(a, b)| *a += *b);
    }
}

#[cfg(test)]
#[test]
fn test_element_wise_add() {
    let mut vec_1: Vec<f64> = vec![0.0, 3.3, 4.4];
    let vec_2: Vec<f64> = vec![1.1, 6.3, -22.4];
    println!("{:?}, {:?}", vec_1, vec_2);
    vec_1.add_assign(&vec_2);
    println!("{:?}, {:?}", vec_1, vec_2);
    let vec_3 = vec_1.add(&vec_2);
    println!("{:?}, {:?}, {:?}", vec_1, vec_2, vec_3);
}
