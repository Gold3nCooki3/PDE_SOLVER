use calculate::*;

/// Performs the Gauss-Seidel / Jacobi algorithm
/// 
/// # Arguments 
/// 
/// * `option` - Defines basic information about the current configuration
/// * `result` - Space to store the results
/// * `arguments` - Stores information about the matrices
/// * `matrix` - three dimensional matrix View where the algorithm is applied on
/// 
/// # Remarks
/// 
/// This is a c style implementation, accessing values with get or uget
/// *Note*: Changing the filed USE_INDEX_CHECKING toggles value checking, with true unsafe code will be executed.
/// 
/// # Panics
/// 
/// If an matrix element is of of bounce, and USE_INDEX_CHECKING ist turn off.
/// 
/// # Errors
/// 
/// If an matrix element is of of bounce, and USE_INDEX_CHECKING ist turn on.
pub fn calc_with_get(options: &Options, results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut Array3<f64>){
    let index_checking : bool = false; 

    // if cfg!(test){ //While Testing index_checking_toggle is copy of mutable basis::USE_INDEX_CHECKING
    //     unsafe{
    //         index_checking = USE_INDEX_CHECKING;
    //     }
    // }else{
    //      index_checking = NO_INDEX_CHECKING;
    // }

    let mut fpisin      : f64 = 0.0;

    let     h                 = arguments.h;
    let mut pih         : f64 = 0.0;
    
    let mut matrix_in   : usize = 0;
    let mut matrix_out  : usize = if options.method == METH_JACOBI{1}else{0};
    let mut term_iteration = options.term_iteration;
    

    if options.inf_func == FUNC_FPISIN
	{
        pih    = PI * h;
		fpisin = 0.25 * TOW_PI_SQARE * h * h;
    }

	while term_iteration > 0
	{
		let mut maxresiduum : f64 = 0.0;

		// over all rows
		for i in 1..arguments.n
		{
			let mut fpisin_i : f64 = 0.0;

			if options.inf_func == FUNC_FPISIN
			{
                fpisin_i = fpisin * f64::sin(pih * i as f64);
            }

			// over all columns
			for j in 1..arguments.n
			{
                let mut star = 
                    if index_checking {
                        0.25 *( matrix[[matrix_in,i-1,j]] + matrix[[matrix_in,i,j-1]] 
                            + matrix[[matrix_in,i,j+1]]+ matrix[[matrix_in,i+1,j]])
                    }else{
                        unsafe{ //ignore index checking
                            0.25 *( matrix.uget((matrix_in,i-1,j)) + matrix.uget((matrix_in,i,j-1))
                                 + matrix.uget((matrix_in,i,j+1)) + matrix.uget((matrix_in,i+1,j)))
                        } 
                };

				if options.inf_func == FUNC_FPISIN {

					star += fpisin_i * f64::sin(pih * j as f64);
				
                }

				if options.termination == TERM_PREC || term_iteration == 1{   
                   
                    let residuum  = 
                        if index_checking {
                            (match matrix.get((matrix_in,i,j)) {
                                                    Some(num) => num,
                                                    None => panic!("Error"),
                                            } - star).abs()
                        }else{
                            unsafe{
                                (matrix.uget((matrix_in,i,j)) - star).abs()
                            }    
                        };

					maxresiduum = maxresiduum.max(residuum);
				}
                
                matrix[[matrix_out,i,j]] = star;
                
			}
        }    
        results.stat_iteration += 1;
        results.stat_precision = maxresiduum;
        
        if options.method == METH_JACOBI {
            let temp = matrix_out;
            matrix_out = matrix_in;
            matrix_in =  temp;
        }
            
		// check for stopping calculation depending on termination method
		if options.termination == TERM_PREC {
			if maxresiduum < options.term_precision {
				term_iteration = 0;
			}
		}
		else if options.termination == TERM_ITER {
			term_iteration -= 1;
		}
    }

    results.m = matrix_in;
}

/// Performs the Gauss-Seidel / Jacobi algorithm
/// 
/// # Arguments 
/// 
/// * `option` - Defines basic information about the current configuration
/// * `result` - Space to store the results
/// * `arguments` - Stores information about the matrices
/// * `matrix_result` - two dimensional matrix View where the final result matrix is copied over
/// 
/// # Remarks
/// 
/// This is a more rust style implementation, accessing the elements with iterators
/// 
/// # Example 
/// 
/// ```text
/// // This does not work since lifetime an mutability differ between matrix_1 matrix_2
/// // even when the lifetime is explicitly set for both to be equal.
/// 
/// fn calc_with_iterators <'a>(options: &CalcOptions, matrix_1 : &'a mut ArrayViewMut2<f64>, matirx_2 : &'a mut ArrayViewMut2<f64>){
///     let star : f64 = 0.0; 
///     let matrix_in  = matrix_1;
///     let matrix_out = matrix_2;
/// 
///     loop{
///         for ((i, j), elt) in matrix_out.indexed_iter_mut(){
///         if i == 0 || j  == 0 || i >= arguments.n || j >= arguments.n {continue}
///         
///              star = 0.25 *( matrix_in[[ i-1, j]]
///                        + matrix_in[[ i, j-1]]
///                        + matrix_in[[ i, j+1]] 
///                        + matrix_in[[ i+1, j]]
///                       );
/// 
///                *elt = star;
///          }
/// 
///          //...
/// 
///          if options.method == METH_JACOBI{
///             mem::swap(&mut matrix_out, &mut matrix_in);
///          //                            ^^^^^^^^^^^^^^
///          //                            lifetime mismatch       
///                                                                
///          }
///     }
/// }
/// ```
pub fn calc_with_iterators(options: &Options, results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut ArrayViewMut3<f64>){
    use std::mem;

    let mut fpisin      : f64 = 0.0;
    let     h                 = arguments.h;
    let mut pih         : f64 = 0.0;
    let mut m           : usize = 1;

    let mut term_iteration = options.term_iteration;

    if options.inf_func == FUNC_FPISIN
	{
        pih    = PI * h;
		fpisin = 0.25 * TOW_PI_SQARE * h * h;
        
    }

    //the split_matrix ensures that matrix_out and matrix_in have the same lifetime
    //this way mem::swap(&mut matrix_out, &mut matrix_in) works
    let (mut matrix_out,mut matrix_in) = matrix.view_mut().split_at(Axis(0),1);
    loop
	{
        if term_iteration <= 0 {break};
		let mut maxresiduum : f64 = 0.0;

		// over all elements
        for ((_, i, j), elt) in matrix_out.indexed_iter_mut(){
            if i == 0 || j  == 0 || i >= arguments.n || j >= arguments.n {continue}
            
            // calculate star 
            let mut star = 0.25 *( matrix_in[[0, i-1, j]]
                         + matrix_in[[0, i, j-1]]
                         + matrix_in[[0, i, j+1]] 
                         + matrix_in[[0, i+1, j]]
                        );


            if options.inf_func == FUNC_FPISIN{

                star += fpisin * f64::sin(pih * i as f64) * f64::sin(pih * j as f64);

            }

            if options.termination == TERM_PREC || term_iteration == 1{   
                
                let residuum = (matrix_in[[0,i,j]] - star).abs();                  

                maxresiduum = maxresiduum.max(residuum);
            }

            if options.method == METH_JACOBI{

                *elt = star;

            }else{

                matrix_in[[0,i,j]] = star;

            }
        }

        results.stat_iteration += 1;
        results.stat_precision = maxresiduum;
        
        if options.method == METH_JACOBI{
            m = if m == 1 {0}else{1};
            mem::swap(&mut matrix_out, &mut matrix_in); //No copying will be performed here
        }
            
		// check for stopping calculation depending on termination method
		if options.termination == TERM_PREC{

			if maxresiduum < options.term_precision{

				term_iteration = 0;
                
			}
		}
		else if options.termination == TERM_ITER{

			term_iteration -= 1;

		}
    }
    results.m = m;
}

/// Calculates the PDE related to the Jacobi algorithm
/// 
/// # Arguments 
/// 
/// * `matrix` - two dimensional matrix view where the algorithm is applied on
/// * `option` - Defines basic information about the current configuration
/// * `result` - Space to store the results
/// * `arguments` - Stores information about the matrices
/// 
/// # Remarks
/// 
/// This is an highly experimental approach using the maximum of high order functions provided by ndarray/ndalgebra
/// 
/// 
///     //  ___ ___ ___ ___ ___    ___ _c_ _c_ _c_ ___      ___ ___ ___ ___ ___ 
///     // |___|_c_|_c_|_c_|___|  |___|_x_|_x_|_x_|___|    |___|___|___|___|___|
///     // |___|_x_|_x_|_x_|___|  |___|_x_|_x_|_x_|___|  c |_x_|_x_|_x_|___|___|  
///     // |___|_x_|_x_|_x_|___|  |___|_x_|_x_|_x_|___|  c |_x_|_x_|_x_|___|___|  
///     // |___|_x_|_x_|_x_|___|  |___|___|___|___|___|  c |_x_|_x_|_x_|___|___|           
///     // |___|___|___|___|___|  |___|___|___|___|___|    |___|___|___|___|___|
///                                                                              
/// The result is calculated by overlapping different matrix views. The four matrix view are defend by the inner values of the matrix,
/// but are shifted in the all directions by one.
/// 
/// # Panics
/// 
/// Panics when method is Gauss-Seidel
/// 
pub fn calc_with_zip(options: &Options, results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut ArrayViewMut2<f64>){
    use ndarray::Zip;
    if options.method == METH_GAUSS_SEIDEL {panic!("This method isn't even close to the Gauß-Seidel algorithm")}
    
    let mut term_iteration = options.term_iteration;
    let     h              = arguments.h;

    //local Arrays
    let mut fpisin_ij   = Array2::<f64>::zeros((arguments.n+1, arguments.n+1));
    let mut scratch     = Array2::<f64>::zeros((arguments.n+1, arguments.n+1));

    if options.inf_func == FUNC_FPISIN{
        let pih    = PI * h;
        let fpisin = 0.25 * TOW_PI_SQARE * h * h;

        //calculates fpisin_ij in advance
        for ((i, j), e)  in fpisin_ij.indexed_iter_mut(){
            *e = fpisin * f64::sin(pih * i as f64) * f64::sin(pih * j as f64)
        }
    }
    
    loop{
        if term_iteration <= 0 {break};
        let mut maxresiduum : f64 = 0.0;
        
        //fill scratch pad
        scratch.fill(0.0);        
        &scratch.slice_mut(s![0.., 0..-1]).zip_mut_with(&matrix.slice(s![0.., 1..]), |y, &n|{*y = *y+n; }); //Shift Matrix left , value from above equal to i - 1
        &scratch.slice_mut(s![1.., 0..]).zip_mut_with(&matrix.slice(s![0..-1, 0..]), |y, &n|{*y = *y+n; }); //Shift Matrix down , value from above equal to j - 1
        &scratch.slice_mut(s![0..-1, 0..]).zip_mut_with(&matrix.slice(s![1.., 0..]), |y, &n|{*y = *y+n; }); //Shift Matrix up , value from above equal to j + 1
        &scratch.slice_mut(s![0.., 1..]).zip_mut_with(&matrix.slice(s![0.., 0..-1]), |y, &n|{*y = *y+n; }); //Shift Matrix right , value from above equal to i + 1

        //transferee relevant values to the matrix
        if options.inf_func == FUNC_F0 {

            matrix.slice_mut(s![1..-1, 1..-1]).zip_mut_with(&scratch.slice(s![1..-1, 1..-1]), |y, &n| {

                let mut residuum = *y;
                
                *y = n * 0.25;

                if options.termination == TERM_PREC || term_iteration == 1{ 

                    residuum = (residuum - *y).abs();

                    maxresiduum = maxresiduum.max(residuum);
                }
            });

        }else{
            //Iterator includes fpisin_ij values
            Zip::from(&mut matrix.slice_mut(s![1..-1, 1..-1]))
                .and(&scratch.slice(s![1..-1, 1..-1]))
                .and(&fpisin_ij.slice(s![1..-1, 1..-1]))
                .apply(|y, &n, &fpi_ij|{

                let mut residuum = *y;
                
                *y = n * 0.25 + fpi_ij;

                if options.termination == TERM_PREC || term_iteration == 1{ 

                    residuum = (residuum - *y).abs();

                    maxresiduum = maxresiduum.max(residuum);
                }
            });

        }
 
        //Set Results
        results.stat_iteration += 1;
        results.stat_precision = maxresiduum;

        // check for stopping calculation depending on termination method
        if options.termination == TERM_PREC{
            if maxresiduum <= options.term_precision{

                term_iteration = -1;
            }
        }else if options.termination == TERM_ITER{

            term_iteration -= 1;
        }
    }
}

pub fn calculate(options: &Options, results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut Array3<f64>){
    let mut fpisin      : f64 = 0.0;
    let     h                 = arguments.h;
    let mut pih         : f64 = 0.0;

    let mut matrix_in   : usize = 0;
    let mut matrix_out  : usize = if options.method == METH_JACOBI{1}else{0};
    let mut term_iteration = options.term_iteration;


    if options.inf_func == FUNC_FPISIN
        {
        pih    = PI * h;
                fpisin = 0.25 * TOW_PI_SQARE * h * h;
    }

        while term_iteration > 0
        {
                let mut maxresiduum : f64 = 0.0;

                // over all rows
                for i in 1..arguments.n
                {
                        let mut fpisin_i : f64 = 0.0;

                        if options.inf_func == FUNC_FPISIN
                        {
                fpisin_i = fpisin * f64::sin(pih * i as f64);
            }

                        // over all columns
                        for j in 1..arguments.n
                        {
                let mut star = unsafe{ //ignore index checking
                            0.25 *( matrix.uget((matrix_in,i-1,j)) + matrix.uget((matrix_in,i,j-1))
                                 + matrix.uget((matrix_in,i,j+1)) + matrix.uget((matrix_in,i+1,j)))
                        };

                                if options.inf_func == FUNC_FPISIN {

                                        star += fpisin_i * f64::sin(pih * j as f64);

                }

                                if options.termination == TERM_PREC || term_iteration == 1{

                    let residuum  = unsafe{
                                (matrix.uget((matrix_in,i,j)) - star).abs()
                            };

                                        maxresiduum = maxresiduum.max(residuum);
                                }

                matrix[[matrix_out,i,j]] = star;

                        }
        }
        results.stat_iteration += 1;
        results.stat_precision = maxresiduum;

        if options.method == METH_JACOBI {
            let temp = matrix_out;
            matrix_out = matrix_in;
            matrix_in =  temp;
        }

                // check for stopping calculation depending on termination method
                if options.termination == TERM_PREC {
                        if maxresiduum < options.term_precision {
                                term_iteration = 0;
                        }
                } else if options.termination == TERM_ITER {
                        term_iteration -= 1;
                }
    }

    results.m = matrix_in;
}