#![allow(unused)]

#[macro_use(s)]
extern crate ndarray;
extern crate ndarray_parallel;
extern crate rayon;
extern crate crossbeam;
extern crate scoped_threadpool;

use ndarray::prelude::*;


pub mod basis; //all structures & constants, methods for in- and output 
    use basis::displaymatrix::display_result;
    use basis::askparams::ask_params;
    use basis::init_matrix::{init_naive};
    use basis::*;

pub mod calculate; //all calculate functions
    use calculate::{serial, parallel};

//just an expensive test function
fn fibonacci_reccursive(n: i32) -> u64 {
	if n < 0 {
		panic!("{} is negative!", n);
	}
	match n {
		0     => panic!("zero is not a right argument to fibonacci_reccursive()!"),
		1 | 2 => 1,
		3     => 2, 
		/*
		50    => 12586269025,
		*/
		_     => fibonacci_reccursive(n - 1) + fibonacci_reccursive(n - 2)
	}
}

fn main() {
    //Define
    use std::time::Instant;
    
    let mut matrix   : Array3<f64>;
    let mut results  = CalcResults::default();
    let mut options  = Options::default(); 
	// by default 1 2 400 1 2 1001
   
    //ask_params(&mut options);  // 3 - 7 sec overhead
    
    let elements = (options.interlines * 8) + 9 - 1;
    let arguments  = CalcArguments {
            n: elements as usize,
            num_matrices: options.method,
            h: 1.0 / elements as f64
    };
    
    {
        matrix = Array::<f64, _>::zeros((2,arguments.n+1, arguments.n+1));
        init_naive(&mut matrix, &options, &arguments);
       
        //Just for debugging ..
        options.number = 1; // <-- number of processes
        options.inf_func = 1; //<-- 2 for a more expensive distubrance function 
        //#################################################
        let mut sec : f64 = 0.0;
        let now = Instant::now(); 

            if options.method == METH_GAUSS_SEIDEL || options.number == 1{
                println!("SERIAL");
                serial::calc_with_get(&options, &mut results, &arguments, &mut matrix); //<-- runtime doubled when calc_get_with_threadpool is in else block 
                //parallel::calc_parallel_light(&options, &mut results, &arguments, &mut  matrix);
               
            }else{
                println!("PARALLEL");
                //serial::calc_with_get(&options, &mut results, &arguments, &mut matrix);
                //parallel::calc_parallel_light(&options, &mut results, &arguments, &mut  matrix);
                parallel::calc_get_with_threadpool(&options, &mut results, &arguments, &mut  matrix);
				//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ if this is comented out the runtime gets reduced to half the orignal runtime .. 				
            }

        let elapsed = now.elapsed();
        sec = (elapsed.as_secs() as f64) + (elapsed.subsec_nanos() as f64 / 1000_000_000.0);      
         
        //#################################################

        results.elapsed_time = sec;
    }
    

display_result(&matrix.subview(Axis(0), results.m), &arguments, &results, &options);
}