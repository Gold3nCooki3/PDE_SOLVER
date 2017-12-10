#![feature(test)]
extern crate pde_solver;
extern crate ndarray;

// #[macro_use] 
extern crate bencher;
extern crate test;

//use bencher::Bencher;
use test::Bencher;

use pde_solver::calculate::parallel;
use pde_solver::basis::{Options, CalcArguments, CalcResults};
use ndarray::prelude::*;

/* ======= jacobi =======*/
#[bench]
fn parallel_jacobi_threadpool(b: &mut Bencher){   
   
    let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 1,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (2,arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

  
    b.iter(|| {parallel::calc_get_with_threadpool(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn parallel_jacobi_crossbeam(b: &mut Bencher){
  
    let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 1,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (2,arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

  
    b.iter(|| {parallel::calc_get_with_crossbeam(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn parallel_jacobi_iter(b: &mut Bencher){
     let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 1,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (2, arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

    b.iter(|| {parallel::calc_with_iterators(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn parallel_jacobi_zip(b: &mut Bencher){
    let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 1,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

    b.iter(|| {parallel::calc_with_zip(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

/* ===================== EXPENSIVE ====================*/

/* ======= jacobi =======*/
#[bench]
#[ignore]
fn parallel_jacobi_threadpool_large(b: &mut Bencher){   
   
    let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 100,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (2,arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

  
    b.iter(|| {parallel::calc_get_with_threadpool(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn parallel_jacobi_crossbeam_large(b: &mut Bencher){
  
    let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 100,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (2,arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

  
    b.iter(|| {parallel::calc_get_with_crossbeam(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn parallel_jacobi_iter_large(b: &mut Bencher){
    let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 100,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (2, arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

    b.iter(|| {parallel::calc_with_iterators(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn parallel_jacobi_zip_large(b: &mut Bencher){
    let mut option      = Options{
            number: 2,
            method: 2,
            interlines: 1000,
            inf_func: 2,
            termination: 2,
            term_precision: 0.0,
            term_iteration: 100,
    };

    let elements = (option.interlines * 8) + 9 - 1;
    let arguments = CalcArguments {
            n: elements as usize,
            num_matrices: option.method + 1,
            h: 1.0 / elements as f64
    };

    let mut results = CalcResults{ ..Default::default()}; 
    let shape = (arguments.n+1, arguments.n+1);
    let mut matrix  = Array::<f64, _>::zeros(shape);

    b.iter(|| {parallel::calc_with_zip(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

// benchmark_group!(small, jacobi_iter, jacobi_wget, jacobi_wuget, jacobi_zip);
// benchmark_group!(benches, jacobi_wuget, jacobi_wget, jacobi_iter, jacobi_zip);
// benchmark_main!(benches);