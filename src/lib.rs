//! This Program calculates the tow dimensional partial differential equations (PDE).
//! This can be done by using the Gauss-Seidel or Jacobi algorithms
//! 
//! The Result will be presented in a simplified view of the matrix and the calculated error coefficient of the equations.
//! 
//! © Nils C. Meier, Universität Hamburg.
//! 
//! # Features 
//! 
//! Compiling with the different features results in different functions for the calculation.
//! 
//! * calc_w_iter  - Accessing matrix elements with an iterator.  
//! * calc_w_zip   - Using a Zip function to add all neighbours and disturbance function in in 5 Operations (only for Jacobi)
//! * default      - Accessing each element using matrix indices
//! 
//! # Usage
//! 
//! To start the program use ./prog [1] [2] [3] [4] [5] [6]
//! 
//! * [1] - Number of process
//! * [2] - Method                :      Gauss-Seidel(1) or Jacobi(2)          
//! * [3] - Interlines            :      (0 .. 1000), will result in (( interlines ⋅ 8 ) + 9 - 1)^2 elements   
//! * [4] - Disturbance function  :      f(x,y) = 0 (1) or f(x,y) = 2π^(2) * sin(π ⋅ x) * sin(π ⋅ y)
//! * [5] - Termination method    :      Sufficient accuracy(1) or Number of iterations(2)
//! * [6] - Error coefficient     :      1e-4 <= x <= 1e-20 
//! * [6] - Number of iterations  :      1    <= x <= 20 000   
//!  
#![allow(unused)]
#[macro_use(s)]

extern crate ndarray;
extern crate itertools;
extern crate ndarray_parallel;
extern crate rayon;
extern crate crossbeam;
extern crate scoped_threadpool;
//extern crate pond;

pub mod basis;
pub mod calculate;