use basis::*;

/// Display statistics an and the final Iteration of the Matrix in an manageable view
///  
/// # Arguments
/// 
/// * `results` - Results to be displayed
/// * `matrix` - Matrix to be displayed 
/// 
/// # Remark 
/// 
/// The Matrix will be shown in an 9 x 9 square view.
/// 
pub fn display_result(matrix: &ArrayView2<f64>, arguments: &CalcArguments, results: &CalcResults, options: &Options){
    display_statistics(arguments, results, options);
    display_matrix(matrix, options);
}

/// Displays the Matrix in an 9 x 9 square
/// 
/// # Arguments 
/// 
/// * `matrix` - An two dimensional view of the result matrix
/// 
/// *Note*: decimal defines the number of decibels shown in the final matrix view.
fn display_matrix(matrix: &ArrayView2<f64>, option: &Options,){
    let     decimal = 8;
    let     pattern : Vec<usize> = (0..9).map(|x| (x * ((*option).interlines + 1) as usize)).collect();
    
    println!("Matrix:");
    let print = matrix.select(Axis(0), &pattern).select(Axis(1), &pattern);
    println!("{:.*}", decimal,print);
}

/// Displays the Configuration and Results of the program run
/// 
/// # Arguments 
/// 
/// * `results` - Containing statistics about the program run
/// * `options` - Containing the configuration of the algorithm 
fn display_statistics(arguments: &CalcArguments, results: &CalcResults, options: &Options){
use std::mem; 
    println!("============================================================
Program for calculation of partial differential equations.
============================================================
(c) Nils C. Meier, Uni Hamburg.\n
============================================================\n");

    let n = arguments.n;

    let method = if options.method == METH_JACOBI {
        "Jacobi"
    }else{
        "Gauss-Seidel"
    };
    let func = if options.inf_func == FUNC_F0 {
        "f(x,y) = 0"
    }else{
        "f(x,y) = 2 π² * sin(π ⋅ x) * sin(π ⋅ y)"
    };
    let term = if options.termination == TERM_PREC{
        "Sufficient accuracy"
    }else{
        "Number of iterations"
    };
    println!("Time for calculation: {} s", results.elapsed_time);
    println!("Used CPU Cores        {}", options.number);
    println!("Memory usage:         {} MiB",  (((2* (n + 1) * (n + 1)) as f64) * ((mem::size_of::<f64>() * arguments.num_matrices as usize) as f64) / 1024.0 / 1024.0));
    println!("Calculating method:   {}", method);
    println!("Interlines:           {}", options.interlines);
    println!("Disturbance function: {}", func);
    println!("Termination method:   {}", term);
    println!("Number of iterations: {}", results.stat_iteration);
    println!("Error coefficient:    {:e}", results.stat_precision);
}
