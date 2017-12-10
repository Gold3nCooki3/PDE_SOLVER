use calculate::*;


/// Calculates the PDE related to the Jacobi algorithm in parallel
///
/// # Arguments
///
/// * `matrix` - Two dimensional matrix View where the algorithm is applied on
/// * `option` - Defines basic information about the current configuration
/// * `result` - Space to store the results
/// * `arguments` - Stores information about the matrices
///
/// # Remarks
///
/// Parallelization is done by the included zip parallel methods in the ndarray_parallel crate.
/// These made it easy to extended the code to run this method in parallel. In order to synchronise maxresiduum
/// a mutex is shared to regulate the access. Since the mutex is borrowed for every element,
/// termination by precision is not realistic to use here and therefor leads to a panic
///
/// To see how the function works look in `calculate::serial::calc_with_zip`
///
/// # Panics
///
/// Panics if the method is Gauss-Seidel or the termination method is termination by precision
pub fn calc_with_zip(options: &Options, results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut ArrayViewMut2<f64>){
    use ndarray::Zip;
    use ndarray_parallel::prelude::*;
    use crossbeam;
    use rayon;
    //use std::sync::{Arc, Mutex};

    if options.method == METH_GAUSS_SEIDEL {panic!("This method is not even close to the Gauß-Seidel algorithm!")}
    if options.termination == TERM_PREC {panic!("This method does not support the termination by precision yet!")}

    let mut term_iteration = options.term_iteration;
    let     h              = arguments.h;

    //local Arrays
    let mut fpisin_ij   = Array2::<f64>::zeros((arguments.n+1, arguments.n+1));
    let mut scratch     = unsafe{Array2::<f64>::uninitialized((arguments.n+1, arguments.n+1))};
    //let maxresiduum_mutex = Arc::new(Mutex::new(0_f64));
    let     maxresiduum_queue : crossbeam::sync::SegQueue<f64>= crossbeam::sync::SegQueue::new() ;
    let     ref_maxresiduum = &maxresiduum_queue;


    //Limit number of threads
    match rayon::initialize(rayon::Configuration::new().num_threads(options.number as usize)){
        Ok(worked) => worked,
        Err(err) => println!("{}", err),
    }

    if options.inf_func == FUNC_FPISIN{
        let pih    = PI * h;
        let fpisin = 0.25 * TOW_PI_SQARE * h * h;

        //per calculating the fpisin_ij
        unsafe{
            let mut array_i   = Array2::<usize>::uninitialized((arguments.n+1, arguments.n+1));
            let mut array_j   = Array2::<usize>::uninitialized((arguments.n+1, arguments.n+1));

            // initialise indices for parallel calculation of 1/4 * 2pi* h^2 * sin(pi * h * i) * sin(pih * j)
            for i in 0 .. arguments.n{
                for j in 0 .. arguments.n{
                    array_i[[i, j]] = i;
                    array_j[[i, j]] = j;
                }
            }

            //calculating parallel
            Zip::from(&mut fpisin_ij)
                .and(&array_i)
                .and(&array_j)
                .par_apply(|e, &i, &j|{
                    *e = fpisin * f64::sin(pih * i as f64) * f64::sin(pih * j as f64)
                });

        } // drop array_i and array_j here
    }

    loop{
        if term_iteration <= 0 {break};

        // {
        //    let mut maxresiduum_shared = *maxresiduum_mutex.lock().unwrap();
        //     maxresiduum_shared = 0.0;
        // } // drop mutex

        //fill scratch pad with 0
        scratch.fill(0.0);

        Zip::from(&mut scratch.slice_mut(s![0.., 0..-1]))
            .and(&matrix.slice(s![0.., 1..]))
            .par_apply(|y, &n|{*y = *y+n; }); //Shift Matrix left , value from above equal to i - 1

        Zip::from(&mut scratch.slice_mut(s![1.., 0..]))
            .and(&matrix.slice(s![0..-1, 0..]))
            .par_apply(|y, &n|{*y = *y+n; }); //Shift Matrix down , value from above equal to j - 1

        Zip::from(&mut scratch.slice_mut(s![0..-1, 0..]))
            .and(&matrix.slice(s![1.., 0..]))
            .par_apply(|y, &n|{*y = *y+n; }); //Shift Matrix up , value from above equal to j + 1

        Zip::from(&mut scratch.slice_mut(s![0.., 1..]))
            .and(&matrix.slice(s![0.., 0..-1]))
            .par_apply(|y, &n|{*y = *y+n; });  //Shift Matrix right , value from above equal to i +

        Zip::from(&mut matrix.slice_mut(s![1..-1, 1..-1]))
            .and(&scratch.slice(s![1..-1, 1..-1]))
            .and(&fpisin_ij.slice(s![1..-1, 1..-1]))
            .par_apply(|y, &n, &fpi_ij|{

            let star =
                if options.inf_func == FUNC_FPISIN {
                    n * 0.25 + fpi_ij
                }else{
                    n * 0.25
                };

            if term_iteration == 1 {
                // accessing the mutex or the parallel Queue for every element is remarkably expensive,
                // so it is only realistic to access it only in the last round
                // TERM_PREC cannot be used here, because it would access the mutex ever round

                let residuum = (*y - star).abs();
                //let mut maxresiduum_shared = maxresiduum_mutex.lock().unwrap();
                //*maxresiduum_shared = (*maxresiduum_shared).max(residuum);
                (*ref_maxresiduum).push(residuum);
            }

            *y = star;
        });

        // access maxresiduum
        //let maxresiduum_shared = maxresiduum_mutex.lock().unwrap();
        //let maxresiduum = *maxresiduum_shared;

        let     maxresiduum : f64 = reduce_result_max(ref_maxresiduum, 0.0);

        // Set Results
        results.stat_iteration += 1;
        results.stat_precision = maxresiduum;

        // Reduce iteration
        term_iteration -= 1;
    }
}

/// Performs the Jacobi algorithm in parallel
///
/// # Arguments
///
/// * `option` - Defines basic information about the current configuration
/// * `result` - Space to store the results
/// * `arguments` - Stores information about the matrices
/// * `matrix_result` - two dimensional matrix view
///
/// # Remarks
///
/// A parallel implementation of the calc_with_get method. This implementation is build on the crossbeam crate to synchronise data. In addition to the code following
/// parts needed to be added
///
/// * `matrix splitting`    - splitting the matrix in p mutable parts to access every part individually
/// * `copying shared rows` - To ensure that no mutable references exits on these rows, it is either necessary to copy all rows every
/// iteration or to split the array in a mutable part and an immutable part while only splitting the mutable part in p parts. This requires to split the array every iteration
/// * `syncing maxresiduum` -
///  Using the JoinHandle structure form the crossbeam create makes it easy to calculate an temporary result in the each process and joining theses afterwards, is not performant
///
/// # Panics
///
/// Panics if Gauss-Seidel method should executed
pub fn calc_get_with_crossbeam(options: &Options, results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut ArrayViewMut3<f64>){
    use std::mem;
    use crossbeam;

    let index_checking : bool = false;

    let mut fpisin      : f64 = 0.0;
    let     h           : f64 = arguments.h;
    let mut pih         : f64 = 0.0;

    let mut matrix_in   : usize = 0;
    let mut matrix_out  : usize = 1;

    let mut term_iteration = options.term_iteration;

    //For parallelization
    let     processes : usize = if options.number as usize <= arguments.n { options.number as usize }else{arguments.n}; // number of processes /currently the main process is idle
    let     size      : usize = arguments.n;
    let     overhead  : usize = (size+1) % processes;

    //initialise Array for all Stripes, shape is [ #processes, 2 rows, #width]
    let mut matrix_stripes = Array::<f64, _>::zeros((processes, 2, arguments.n+1));

    //initialise parallel Queue for sync of maxresiduum
    let     maxresiduum_queue : crossbeam::sync::SegQueue<f64>= crossbeam::sync::SegQueue::new() ;
    let     ref_maxresiduum = &maxresiduum_queue;  //reference can be copied

    //first split of the array, initialise of vector of all mutable views
    let (left , right) = matrix.view_mut().split_at(Axis(1),(size+1)/processes);
    let mut matrix_blocks = vec![left];
    //splitting the array in more blocks by splitting the right part again,
    {
        let mut temp = right;
        for i in 1 .. processes{
            let spilt_ix = if i < processes-1 {(size+1)/processes} else{(size+1)/processes + overhead} ;
            let (left , new_right) = temp.split_at(Axis(1),spilt_ix);
            matrix_blocks.push(left);
            temp = new_right // last right block is empty, because the split is on the last row ( + overhead)
        }
    } //temp is dropped here

    if options.inf_func == FUNC_FPISIN
	{
        pih    = PI * h;
		fpisin = 0.25 * TOW_PI_SQARE * h * h;
    }

    while term_iteration > 0 {
        let mut maxresiduum : f64 = 0.0;

        copy_matrix_stripes(&matrix_blocks, matrix_in, processes, &mut matrix_stripes.view_mut());
        let ref_matrix_stripes = &matrix_stripes;

        crossbeam::scope(|scope| {
            for (__thread_id, matrix) in matrix_blocks.iter_mut().enumerate(){
                scope.spawn(move || {

                    //set up start and end
                    let start = if __thread_id > 0 {0}else{1}; // first process will skip first row
                    let end = if __thread_id < processes-1 {((size+1)/processes)}else{((size+1)/processes)+overhead-1}; //last process takes overhead

                    // over all rows
                    for i in start..end
                    {
                        let mut fpisin_i : f64 = 0.0;

                        if options.inf_func == FUNC_FPISIN
                        {
                            fpisin_i = fpisin * f64::sin(pih * (i + __thread_id * (size+1)/processes) as f64);
                        }

                        // over all columns
                        for j in 1..arguments.n
                        {
                            let mut star : f64 =
                                if index_checking {
                                    0.25 *(
                                    if i > 0 { matrix[[matrix_in, i-1,j]] }
                                        else{ ref_matrix_stripes[[__thread_id-1, 1, j]] } // uses the last row of the process -1 (above in the matrix)
                                    + matrix[[matrix_in, i,j-1]]
                                    + matrix[[matrix_in, i,j+1]]
                                    + if i < end - 1 || __thread_id == processes-1  { matrix[[matrix_in, i+1,j]] }
                                        else{ ref_matrix_stripes[[__thread_id+1, 0, j]]})  // uses the first row of the process + 1 (below in the matrix)
                                }else{
                                    unsafe{
                                        0.25 *(
                                            if i > 0 {
                                                    *matrix.uget((matrix_in,i-1,j))
                                                }else{
                                                    ref_matrix_stripes[[__thread_id-1, 1, j]]
                                                }
                                            + matrix.uget((matrix_in,i,j-1))
                                            + matrix.uget((matrix_in,i,j+1))
                                            + if i < end - 1 || __thread_id == processes-1  {
                                                    *matrix.uget((matrix_in,i+1,j))
                                                }else{
                                                    ref_matrix_stripes[[__thread_id+1, 0, j]]
                                            })
                                    }
                                };

                            if options.inf_func == FUNC_FPISIN{
                                star += fpisin_i * f64::sin(pih * j as f64);
                            }

                            if options.termination == TERM_PREC || term_iteration == 1{
                                let residuum : f64 = if index_checking {
                                    (matrix[[matrix_in,i,j]] - star).abs()
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

                    (*ref_maxresiduum).push(maxresiduum);
                });
            }
        });

        maxresiduum = reduce_result_max(ref_maxresiduum, 0.0);

        results.stat_iteration += 1;
        results.stat_precision = maxresiduum;

        if options.method == METH_JACOBI{
            mem::swap(&mut matrix_in, &mut matrix_out);
        }else{
            panic!{"The Gaus Seidel Algorthm cannot be calculatet like this!"}
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

/// Performs the Jacobi algorithm in parallel
///
/// # Arguments
///
/// * `option` - Defines basic information about the current configuration
/// * `result` - Space to store the results
/// * `arguments` - Stores information about the matrices
/// * `matrix_result` - two dimensional matrix view.
///
/// # Remarks
///
/// A parallel implementation of the calc_with_get method. This implementation is build on the thread pool crate to reduce the overhead
/// with crating new threads in very iteration. In addition to the code following parts needed to be added
///
/// * `matrix splitting`    - splitting the matrix in p mutable parts to access every part individually
/// * `copying shared rows` - To ensure that no mutable references exits on these rows, it is either necessary to copy all rows every
/// iteration or to split the array in a mutable part and an immutable part while only splitting the mutable part in p parts. This requires to split the array every iteration
/// * `syncing maxresiduum` - Mutex to sync the maxresiduum at the end of each parallel block
///
/// # Panics
///
/// Panics if Gauss-Seidel method should executed
///
/// # Safety
///
/// Due to a bug in the crate this version might be leaking memory https://github.com/Kimundi/scoped-threadpool-rs/issues/16
/// As an alternative pond was tested. The memory leak count be proofed, but pond seemed slower than scoped-threadpool crate
pub fn calc_get_with_threadpool(options: &Options,  results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut Array3<f64>){
    use std::mem;
    use std::sync::{Arc, Mutex};
    use scoped_threadpool;
    //use pond;
    use std::time::Instant;

    let mut fpisin      : f64 = 0.0;
    let     h           : f64 = arguments.h;
    let mut pih         : f64 = 0.0;
    let mut matrix_in   : usize = 0;
    let mut matrix_out  : usize = 1;
    let mut term_iteration = options.term_iteration;

    //For parallelization
    let     processes : usize = if options.number as usize <= arguments.n { options.number as usize }else{arguments.n}; // number of processes /currently the main process is idle
    let mut pool              = scoped_threadpool::Pool::new(processes as u32);
    //let mut pool              = pond::Pool::new(processes)
    let     size      : usize = arguments.n;
    let     overhead  : usize = (size+1) % processes;

    //initialise mutex for synchronising maxresiduum between the processes
    let     maxresiduum_mutex = Arc::new(Mutex::new(0_f64));
    let     ref_maxresiduum_mutex = &maxresiduum_mutex;

    //initialise array for all Stripes, shape is [ #processes, 2 rows, #width]
    let mut matrix_stripes = Array::<f64, _>::zeros((processes, 2, arguments.n+1));

    //first split of the array, initialise of vector of all mutable views
    let (left , right) = matrix.view_mut().split_at(Axis(1),(size+1)/processes);
    let mut matrix_blocks = vec![left];
    //splitting the array in more blocks by splitting the right part again,
    {

        let mut temp = right;
        for i in 1 .. processes{

            let spilt_ix = if i < processes-1 {(size+1)/processes} else{(size+1)/processes + overhead} ;
            let (left , new_right) = temp.split_at(Axis(1),spilt_ix);
            matrix_blocks.push(left);
            temp = new_right // last right block is empty, because the split is on the last row ( + overhead)
        }
    } 

    if options.inf_func == FUNC_FPISIN {

        pih    = PI * h;
		fpisin = 0.25 * TOW_PI_SQARE * h * h;
    }

    while term_iteration > 0 {
        { // block for mutex
            let mut shared_maxresiduum = (*ref_maxresiduum_mutex).lock().unwrap();
            *shared_maxresiduum = 0.0;
            // mutex goes out of scope an will be unlocked
        }
        //copy the shared rows between to processes into matrix_stripes, to safely the rows
        copy_matrix_stripes(&matrix_blocks, matrix_in, processes, &mut matrix_stripes.view_mut());
        let ref_matrix_stripes = &matrix_stripes; //bind ownership to immutable reference to accessed the matrix in the parallel block

        pool.scoped( |scoped| {
            for (__thread_id, block) in matrix_blocks.iter_mut().enumerate(){
                scoped.execute(move || {

                    let mut maxresiduum : f64 = 0.0;
                    //set up start and end
                    let start = if __thread_id > 0 {0}else{1}; // first process will skip first row
                    let end = if __thread_id < processes-1 {((size+1)/processes)}else{((size+1)/processes)+overhead-1}; //last process takes overhead

                    for i in start .. end{

                        let mut fpisin_i : f64 = 0.0;

                        if options.inf_func == FUNC_FPISIN {
                            fpisin_i = fpisin * f64::sin(pih * (i + __thread_id * (size+1)/processes) as f64);
                        }

                        for j in 1 .. size{
                            let mut star =  unsafe { 0.25 *(
                                if i > 0 { *block.uget((matrix_in, i-1,j)) }
                                    else{ 0.0 }//ref_matrix_stripes[[__thread_id-1, 1, j]] } // uses the last row of the process -1 (above in the matrix)
                                + *block.uget((matrix_in, i,j-1))
                                + *block.uget((matrix_in, i,j+1))
                                + if i < end - 1 || __thread_id == processes-1  { *block.uget((matrix_in, i+1,j)) }
                                    else{ 0.0 })// ref_matrix_stripes[[__thread_id+1, 0, j]]}) ; // uses the first row of the process + 1 (below in the matrix)
                            };
                            if options.inf_func == FUNC_FPISIN {
                                star += fpisin_i * f64::sin(pih * j as f64);
                            }

                            if options.termination == TERM_PREC || term_iteration == 1 {
                                let residuum = (block[[matrix_in,i,j]] - star).abs();
                                maxresiduum = maxresiduum.max(residuum);
                            }

                            block[[matrix_out,i,j]] = star;
                        }
                    }
                    let mut shared_maxresiduum = (*ref_maxresiduum_mutex).lock().unwrap();  // lock mutex
                    *shared_maxresiduum = maxresiduum.max(*shared_maxresiduum);             // write maxresiduum
                });
            }
        });

        let shared_maxresiduum = (*ref_maxresiduum_mutex).lock().unwrap();

        results.stat_iteration += 1;
        results.stat_precision = *shared_maxresiduum;

        if options.method == METH_JACOBI{
            mem::swap(&mut matrix_out, &mut matrix_in);
        }else{
            panic!{"The Gaus Seidel Algorthm cannot be calculatet like this!"}
        }

        // check for stopping calculation depending on termination method
        if options.termination == TERM_PREC {
            if *shared_maxresiduum < options.term_precision {
                term_iteration = 0;
            }
        }
        else if options.termination == TERM_ITER {
            term_iteration -= 1;
        }
    }
    drop(pool);
    results.m = matrix_in;
}


/// Performs the Jacobi algorithm in parallel
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
/// Replacing the iterator with an Zip + par_apply and an mutex for synchronising the maxresiduum made the it possible to execute the code in parallel.
/// The indexed iterator does not support the functionality to be run in parallel in this version of the creates.
///
/// # Panics
///
/// Panics if the method is Gauss-Seidel or the termination method is termination by precision
pub fn calc_with_iterators(options: &Options, results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut ArrayViewMut3<f64>){
    use std::mem;
    use ndarray::Zip;
    use ndarray_parallel::prelude::*;
    //use std::sync::{Arc,Mutex};
    use crossbeam;
    use rayon;

    if options.termination == TERM_PREC {panic!("This method does not support the termination by precision yet!")}

    let mut fpisin      : f64 = 0.0;

    let     h                 = arguments.h;
    let mut pih         : f64 = 0.0;
    let mut m           : usize = 1;

    let     matrix_in_ix: usize = 0;
    let mut term_iteration = options.term_iteration;

    //let     maxresiduum_mutex = Arc::new(Mutex::new(0_f64));
    let     maxresiduum_queue : crossbeam::sync::SegQueue<f64>= crossbeam::sync::SegQueue::new() ;
    let     ref_maxresiduum = &maxresiduum_queue;

    //Limit number of threads
    match rayon::initialize(rayon::Configuration::new().num_threads(options.number as usize)){
        Ok(worked) => worked,
        Err(err) => println!("{}", err),
    }

    //make array only for indexing
    let mut array_i   = unsafe{ Array2::<usize>::uninitialized((arguments.n+1, arguments.n+1))};
    let mut array_j   = unsafe{ Array2::<usize>::uninitialized((arguments.n+1, arguments.n+1))};

    if options.inf_func == FUNC_FPISIN
	{
        pih    = PI * h;
		fpisin = 0.25 * TOW_PI_SQARE * h * h;

    }

    for i in 0 .. arguments.n{
        for j in 0 .. arguments.n{
            array_i[[i, j]] = i;
            array_j[[i, j]] = j;
        }
    }

    //the split_matrix ensures that matrix_out and matrix_in have the same lifetime
    //this way mem::swap(&mut matrix_out, &mut matrix_in) works
    let (mut matrix_out,mut matrix_in) = matrix.view_mut().split_at(Axis(0),1);
    loop
	{
        if term_iteration <= 0 {break};

        //{ // block for mutex
        //    let mut shared_maxresiduum = maxresiduum_mutex.lock().unwrap();
        //    *shared_maxresiduum = 0.0;
        //    // mutex goes out of scope an will be unlocked
        //}

		// over all elements
        Zip::from(&mut matrix_out.subview_mut(Axis(0), 0))
            .and(&array_i)
            .and(&array_j)
            .par_apply(| elt, &i, &j| {
                if i == 0 || j  == 0 || i >= arguments.n || j >= arguments.n { }else{

                    // calculate star
                    let mut star = 0.25 *( matrix_in[[matrix_in_ix, i-1, j]]
                                + matrix_in[[matrix_in_ix, i, j-1]]
                                + matrix_in[[matrix_in_ix, i, j+1]]
                                + matrix_in[[matrix_in_ix, i+1, j]]
                                );


                    if options.inf_func == FUNC_FPISIN{

                        star += fpisin * f64::sin(pih * i as f64) * f64::sin(pih * j as f64);
                    }

                    if term_iteration == 1{
                        // accessing the mutex for every element is remarkably expensive,
                        // so it is only realistic to access it only in the last round
                        // TERM_PREC cannot be used here, because it would access the mutex ever round

                        let residuum = (matrix_in[[matrix_in_ix,i,j]] - star).abs();
                        //let mut shared_maxresiduum = maxresiduum_mutex.lock().unwrap();  // lock mutex
                        //*shared_maxresiduum = shared_maxresiduum.max(residuum);           // write maxresiduum
                        (*ref_maxresiduum).push(residuum);
                    }

                    if options.method == METH_JACOBI{

                        *elt = star;
                    }else{
                        panic!("The Gaus Seidel Algoritm cannot be performed here");
                        //matrix_in[[matrix_in_ix,i,j]] = star;
                        //^^^^^^^^^ cannot be borrowed mutably
                    }
                }
            });

        //let shared_maxresiduum = maxresiduum_mutex.lock().unwrap();  // lock mutex
        //maxresiduum = *shared_maxresiduum;

        let     maxresiduum : f64 = reduce_result_max(ref_maxresiduum, 0.0);

        results.stat_precision = maxresiduum;
        results.stat_iteration += 1;

        if options.method == METH_JACOBI{
            m = if m == 1 {0}else{1};
            mem::swap(&mut matrix_out, &mut matrix_in); //No copying will be performed here
        }

		// check for stopping calculation depending on termination method
		term_iteration -= 1;
    }
    results.m = m;
}


pub fn calc_parallel_light(options: &Options,  results: &mut CalcResults, arguments: &CalcArguments, matrix: &mut Array3<f64>){
    use std::mem;
    use ndarray_parallel::prelude::*;
    use std::sync::{Arc, Mutex};
    use scoped_threadpool;
    use rayon;
    use ndarray::Zip;    //use pond;
    use std::time::Instant;

    let mut fpisin      : f64 = 0.0;
    let     h           : f64 = arguments.h;
    let mut pih         : f64 = 0.0;
    let mut matrix_in   : usize = 0;
    let mut matrix_out  : usize = 1;
    let mut term_iteration = options.term_iteration;
    

    //For parallelization
    let     processes : usize = if options.number as usize <= arguments.n { options.number as usize }else{arguments.n}; // number of processes /currently the main process is idle
    let mut pool              = scoped_threadpool::Pool::new(processes as u32);
    //let mut pool              = pond::Pool::new(processes)
    let     size      : usize = arguments.n;
    let     overhead  : usize = (size+1) % processes;
    let part  = (size+1)/processes;
    let part_oh = part + overhead-1;

    //initialise mutex for synchronising maxresiduum between the processes
    let     maxresiduum_mutex = Arc::new(Mutex::new(0_f64));
    let     ref_maxresiduum_mutex = &maxresiduum_mutex;
    // let     maxtrix_mutex = Arc::new(Mutex::new(matrix));
    // let     ref_maxtrix_mutex = &maxtrix_mutex;
    // let now = Instant::now();

    // first split of the array, initialise of vector of all mutable views
    let (left , right) = matrix.view_mut().split_at(Axis(1),(size+1)/processes);
    
    let mut matrix_blocks = vec![left];
    //splitting the array in more blocks by splitting the right part again,
    if processes > 1{

        let mut temp = right;
        for i in 1 .. processes{

            let spilt_ix = if i < processes-1 {(size+1)/processes} else{(size+1)/processes + overhead} ;
            let (left , new_right) = temp.split_at(Axis(1),spilt_ix);
            matrix_blocks.push(left);
            temp = new_right // last right block is empty, because the split is on the last row ( + overhead)
        }
    } //temp is dropped here
     
    // let elapsed = now.elapsed();
    // let sec = (elapsed.as_secs() as f64) + (elapsed.subsec_nanos() as f64 / 1000_000_000.0);
    // println!("{}", sec);

    let mut fpisin_ij   = unsafe{Array2::<f64>::uninitialized((size+1, size+1))};
    

    if options.inf_func == FUNC_FPISIN {

        match rayon::initialize(rayon::Configuration::new().num_threads(options.number as usize)){
            Ok(worked) => worked,
            Err(_) => (),
        }

        pih    = PI * h;
		fpisin = 0.25 * TOW_PI_SQARE * h * h;

       unsafe{
            let mut array_i   = Array2::<usize>::uninitialized((arguments.n+1, arguments.n+1));
            let mut array_j   = Array2::<usize>::uninitialized((arguments.n+1, arguments.n+1));

            // initialise indices for parallel calculation of 1/4 * 2pi* h^2 * sin(pi * h * i) * sin(pih * j)
            for i in 0 .. arguments.n{
                for j in 0 .. arguments.n{
                    array_i[[i, j]] = i;
                    array_j[[i, j]] = j;
                }
            }

            //calculating parallel
            Zip::from(&mut fpisin_ij)
                .and(&array_i)
                .and(&array_j)
                .par_apply(|e, &i, &j|{
                    *e = fpisin * f64::sin(pih * i as f64) * f64::sin(pih * j as f64)
                });

        } 
    }
    let fpisin_ij_ref = &fpisin_ij;

    while term_iteration > 0 {
        if options.termination == TERM_PREC || term_iteration == 1 { // block for mutex
            let mut shared_maxresiduum = (*ref_maxresiduum_mutex).lock().unwrap();
            *shared_maxresiduum = 0.0;
            // mutex goes out of scope an will be unlocked
        }

        pool.scoped( |scoped| {
            for (__thread_id, block) in matrix_blocks.iter_mut().enumerate(){
                scoped.execute(move || {
                    //let mut block = (*ref_maxtrix_mutex).lock().unwrap();
                    let mut maxresiduum : f64 = 0.0;
                    //set up start and end
                    let start = if __thread_id > 0 {0}else{1}; // first process will skip first row
                    let end = if __thread_id < processes-1 {part}else{part_oh}; //last process takes overhead6
                    for i in start .. end{

                        for j in 1 .. size{
                            let mut star =  unsafe { 0.25 *(
                                if i > 0 { *(*block).uget((matrix_in, i-1,j)) }
                                    else{ 0.0 }
                                + *(*block).uget((matrix_in, i,j-1))
                                + *(*block).uget((matrix_in, i,j+1))
                                + if i < end - 1 || __thread_id == processes-1  { *(*block).uget((matrix_in, i+1,j)) }
                                    else{ 0.0 })
                            };

                            if options.inf_func == FUNC_FPISIN {
                                star += unsafe{ fpisin_ij_ref.uget((i + __thread_id * end, j)) };
                            }

                            if options.termination == TERM_PREC || term_iteration == 1 {
                                let residuum = (block[[matrix_in,i,j]] - star).abs();
                                maxresiduum = maxresiduum.max(residuum);

                                let mut shared_maxresiduum = (*ref_maxresiduum_mutex).lock().unwrap();  // lock mutex
                                *shared_maxresiduum = maxresiduum.max(*shared_maxresiduum);             // write maxresiduum
                            }

                            unsafe{
                                *block.uget_mut((matrix_out,i,j))= star;
                            }
                        }
                    }
                });
            }
        });

        let shared_maxresiduum = if options.termination == TERM_PREC || term_iteration == 1 {
             *(*ref_maxresiduum_mutex).lock().unwrap()
        }else{
            0.0
        };

        results.stat_iteration += 1;
        results.stat_precision = shared_maxresiduum;

        if options.method == METH_JACOBI{
            mem::swap(&mut matrix_out, &mut matrix_in);
        }else{
            panic!{"The Gaus Seidel Algorthm cannot be calculatet like this!"}
        }

        // check for stopping calculation depending on termination method
        if options.termination == TERM_PREC {
            if shared_maxresiduum < options.term_precision {
                term_iteration = 0;
            }
        }
        else if options.termination == TERM_ITER {
            term_iteration -= 1;
        }
    }
    drop(pool);
    results.m = matrix_in;
}