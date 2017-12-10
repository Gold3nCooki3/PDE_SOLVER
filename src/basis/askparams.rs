use std;
use basis::*;

/// Validates a given value dependent on the index.
/// 
/// # Arguments 
/// 
/// * `arg` - Argument to check
/// * `index` - defines the range of the value
/// 
/// # Panics
/// When the index is wrong.
fn test_params( arg:&i32, index:usize ) -> bool{
    let argv = *arg;
    match index{
        0 => argv > 0,
        1 => argv == METH_GAUSS_SEIDEL || argv == METH_JACOBI,
        2 => argv != 1 && argv <= 10000,
        3 => argv == FUNC_F0    || argv == FUNC_FPISIN,
        4 => argv == TERM_PREC  || argv == TERM_ITER,
        5 => argv > 0 && argv < MAX_ITERATION, 
        _ => panic!("Unexpected Outcome")
    }
}

/// Collects all command line arguments an matches them into a given `Options` structure
/// 
/// # Arguments 
/// 
/// * `option` - Reference to option structure where the all arguments are stored.
/// 
/// # Panics
/// 
/// When the in put string does not match the the requirements.    
pub fn ask_params(option: &mut Options){
    let args: Vec<String>    = std::env::args().collect();

    match args.len(){
    7 =>{ 
            let opt : Vec<i32> = args.iter().skip(1).enumerate().filter(|&(i,_)| i < 5).map(|(i, x)| {
                match x.parse::<i32>(){
                    Ok(num) => {
                        if test_params(&num, i){
                            num
                        }else{
                            panic!("Error: Please fill out correctly!");
                        }
                    },
                    Err(error) => panic!("Error {}", error),
                }
            },)
            .collect();

            if opt[4] == TERM_PREC {
                (*option).term_iteration = MAX_ITERATION;
                (*option).term_precision = match args[6].parse::<f64>(){
                    Ok(num) => if num > 1e-20 && num < 1e-4 
                    {
                       num
                    }else{
                        panic!("Error: Please fill out correctly!");
                    },
                    Err(error) => panic!("Error {}", error),
                }
            }else{
               (*option).term_precision = 0.0;
               (*option).term_iteration = match args[6].parse::<i32>(){
                    Ok(num) => if test_params(&num, 5) 
                    {
                        num
                    }else{
                        panic!("Error: Please fill out correctly!");
                    },
                    Err(error) => panic!("Error {}", error),
                }
            };
            
            (*option).number      = opt[0];
            (*option).method      = opt[1];
            (*option).interlines  = opt[2] as u32;
            (*option).inf_func  = opt[3];
            (*option).termination = opt[4];
        }
    _ => panic!("PANIC"),
    }
}