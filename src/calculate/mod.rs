//! The calculate modules provides methods to solve tow dimensional partial differential equations.
//! 
//! # Serial
//! 
//! Using the Gauss-Seidel or Jacobi algorithm to compute the PDE 
//! 
//! # Parallel
//! 
//! Implementations of the Jacobi algorithm running in parallel

use basis::*;
use ndarray::prelude::*;
use ndarray::Axis;
use crossbeam::sync;

pub mod serial;
pub mod parallel;

#[doc(hidden)]
static mut USE_INDEX_CHECKING : bool = true;

#[doc(hidden)]
pub fn change_index_checking(val: bool){
    unsafe{
        USE_INDEX_CHECKING = val;
    }
}

/// copying values - references are not valid here, the other process has a mutable reference over this part respectively
/// copy the shared rows between to processes into matrix_stripes, to safely the rows
/// 
#[allow(unused)]
fn copy_matrix_stripes<'a, 'b>(matrix_blocks: &'b Vec<ArrayViewMut3<f64>>, matrix_in : usize, processes: usize, matrix_stripes: &'a mut ArrayViewMut3<f64>){

    for p in 0 .. processes{
        // get view of first and last row for each array part
        let last_row_ix = matrix_blocks[p].shape()[1] - 1; // last row is defined by rows - 1 | shape returns &[matrix, rows, columns] 
        let stripe = matrix_blocks[p].subview(Axis(0), matrix_in).select(Axis(0), &[0, last_row_ix]); // select first and last row
        for ((i, j), e) in stripe.indexed_iter(){
            matrix_stripes[[p, i, j]] = *e; // copy each row in the stripe matrix
        }
    }
}

///
/// 
/// 
#[allow(unused)]
fn reduce_result_max(ref_maxresiduum: &sync::SegQueue<f64>, default: f64) -> f64{
    let mut result : f64 = default;
    loop {
        match (*ref_maxresiduum).try_pop() {
            Some(val) => result = result.max(val),
            None => break,
        }
    }
    result
}

///
/// 
/// 
#[allow(unused)]
fn reduce_result_min(ref_maxresiduum: &sync::SegQueue<f64>, default: f64) -> f64{
    let mut result : f64 = default;
    loop {
        match (*ref_maxresiduum).try_pop() {
            Some(val) => result = result.min(val),
            None => break,
        }
    }
    result
}


///
/// 
/// 
#[doc(hidden)]
#[allow(unused)]
fn expensive_fn(a: usize, b: usize) -> usize{
    //calculate Greatest Common Factor
    if b == 0 {
        a
    }else if a == 0{
        b
    }else if a > b{
        expensive_fn(a-b, b)
    }else{
        expensive_fn(a, b-a)
    }
}

