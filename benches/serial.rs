#![feature(test)]
extern crate pde_solver;
extern crate ndarray;

// #[macro_use] 
extern crate bencher;
extern crate test;

//use bencher::Bencher;
use test::Bencher;

use pde_solver::calculate::{serial,change_index_checking};
use pde_solver::basis::{Options, CalcArguments, CalcResults};
use ndarray::prelude::*;

/* ======= gaus_seidel =======*/
#[bench]
fn serial_gaus_seidel_wget(b: &mut Bencher){
    change_index_checking(true);
    
    let mut option      = Options{
            number: 1,
            method: 1,
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

    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn serial_gaus_seidel_uget(b: &mut Bencher){
    change_index_checking(false);
    
    let mut option      = Options{
            number: 1,
            method: 1,
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

    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn serial_gaus_seidel_iter(b: &mut Bencher){    
    let mut option      = Options{
            number: 1,
            method: 1,
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

    b.iter(|| {serial::calc_with_iterators(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

/* ======= jacobi =======*/
#[bench]
fn serial_jacobi_wget(b: &mut Bencher){
    change_index_checking(true);    
   
    let mut option      = Options{
            number: 1,
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

  
    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn serial_jacobi_wuget(b: &mut Bencher){
    change_index_checking(false);
  
    let mut option      = Options{
            number: 1,
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

  
    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn serial_jacobi_iter(b: &mut Bencher){
     let mut option      = Options{
            number: 1,
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

    b.iter(|| {serial::calc_with_iterators(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
fn serial_jacobi_zip(b: &mut Bencher){
     let mut option      = Options{
            number: 1,
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

    b.iter(|| {serial::calc_with_zip(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

/* ===================== EXPENSIVE ====================*/

/* ======= jacobi =======*/
#[bench]
#[ignore]
fn serial_serial_jacobi_wget_large(b: &mut Bencher){
    change_index_checking(true);
   
    let mut option      = Options{
            number: 1,
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

  
    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn serial_jacobi_wuget_large(b: &mut Bencher){
    change_index_checking(false);
    let mut option      = Options{
            number: 1,
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

    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn serial_jacobi_iter_large(b: &mut Bencher){
    let mut option      = Options{
            number: 1,
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

    b.iter(|| {serial::calc_with_iterators(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn serial_jacobi_zip_large(b: &mut Bencher){
    let mut option      = Options{
            number: 1,
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

    b.iter(|| {serial::calc_with_zip(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

/* ======= gaus_seidel =======*/
#[bench]
#[ignore]
fn serial_gaus_seidel_wget_large(b: &mut Bencher){
    change_index_checking(true);    
    
    let mut option      = Options{
            number: 1,
            method: 1,
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

  
    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn serial_gaus_seidel_uget_large(b: &mut Bencher){
    change_index_checking(false);
    
    let mut option      = Options{
            number: 1,
            method: 1,
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

  
    b.iter(|| {serial::calc_with_get(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}

#[bench]
#[ignore]
fn serial_gaus_seidel_iter_large(b: &mut Bencher){
    
    let mut option      = Options{
            number: 1,
            method: 1,
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
 
    b.iter(|| {serial::calc_with_iterators(&mut option, &mut results, &arguments, &mut matrix.view_mut())});
}



// benchmark_group!(gaus_seidel,   gaus_seidel_wget_large, gaus_seidel_uget_large, gaus_seidel_iter_large);
// benchmark_group!(jacobi,        jacobi_wuget_large,     jacobi_wget_large,      jacobi_iter_large,      jacobi_zip_large);
// benchmark_group!(small, jacobi_iter, jacobi_wget, jacobi_wuget, jacobi_zip);
// benchmark_group!(expensive,     jacobi_wuget_large,     jacobi_wget_large,      gaus_seidel_uget_large, gaus_seidel_wget_large, jacobi_zip_large);
// benchmark_group!(benches,       gaus_seidel_wget,       gaus_seidel_uget,       gaus_seidel_iter,       jacobi_wuget,           jacobi_wget,    jacobi_iter,    jacobi_zip);
// benchmark_main!(benches);