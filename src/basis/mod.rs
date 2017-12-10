//! The Basic module contains:
//! 
//! * Constants  
//! * Data structures containing the state of the program
//! * Display Matrix functions
//! * Import of the command line arguments
use ndarray::prelude::*;
use ndarray::Axis;
use std::f64::consts;


pub const MAX_ITERATION     : i32 = 200000;
pub const METH_GAUSS_SEIDEL : i32 = 1;
pub const METH_JACOBI 		: i32 = 2;
pub const FUNC_F0			: i32 = 1;
pub const FUNC_FPISIN		: i32 = 2;
pub const TERM_PREC		    : i32 = 1;
pub const TERM_ITER         : i32 = 2;
pub const TOW_PI_SQARE      : f64 = consts::PI *  consts::PI * 2.0;
pub const PI                : f64 = consts::PI;

/// Contains information concerning the algorithm
/// Can be initialized by the Default function
/// # Fields
/// 
/// * `number`  - number of processes running the algorithm
/// * `method`  - Gauss-Seidel(1) or Jacobi(2)
/// * `interlines`  - number of interlines -> *9+8 lines in the matrix
/// * `inf_func`  - f(x,y) = 0 (1) or f(x,y) = 2pi^2*sin(pi*x)sin(pi*y) (2)
/// * `termination`  - Termination by precision (1) or after a set amount of interlines (2)
/// * `term_iteration`  - interlines till termination
/// * `term_precision`  - precision limit for termination 
#[derive(Debug)]
#[derive(Default)]
pub struct Options {
    pub number          : i32,   /* Number of threads                              */
    pub method          : i32,   /* Gauss Seidel or Jacobi method of iteration     */
    pub interlines      : u32,   /* matrix size = interlines*8+9                   */
    pub inf_func        : i32,   /* inference function                             */
    pub termination     : i32,   /* termination condition                          */
    pub term_iteration  : i32,   /* terminate if iteration number reached          */
    pub term_precision  : f64,   /* terminate if precision reached                 */
}  
    impl Options {
        /// Initialises an Option structure with 0 in all fields
        pub fn default() -> Options {
            Options { number:  1,
                    method: 2,
                    interlines: 400,
                    inf_func: 1,
                    termination: 2,
                    term_iteration: 1001,
                    term_precision: 0.0,
            }
        }
    }

/// Contains information concerning the matrices
/// 
/// # Fields
/// 
/// * `n` - number of lines in the matrix (in both dimensions)
/// * `num_matrices` - number of matrices used dependent on Method 
/// * `h` - change rate by initialisation 
#[derive(Debug)]
pub struct CalcArguments {
    pub n               : usize,
    pub num_matrices    : i32,
    pub h               : f64,        
}

/// Contains the result
/// 
/// # Fields 
/// 
/// * `m` - matrix where the the last iteration is stored
/// * `stat_iteration` - performed iterations 
/// * `stat_precision` - maximum error
#[derive(Debug)]
#[derive(Default)]
pub struct CalcResults {
    pub m               : usize,
    pub stat_iteration  : i32,
    pub stat_precision  : f64,
    pub elapsed_time    : f64,        
}
    impl CalcResults {
        pub fn default() -> CalcResults {
            CalcResults { m: 0,
                    stat_iteration: 0,
                    stat_precision: 0.0,
                    elapsed_time: 0.0,
            }
        }
    }

pub mod displaymatrix;
pub mod askparams;
pub mod init_matrix;