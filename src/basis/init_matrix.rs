use basis::*;


/// Fills one stripe of elements with an constantly dropping / rising value 
/// 
/// # Arguments 
///    
/// * `matrix` - A two dimensional Array View that will be initialized
/// * `h`      - Change rate of value
/// * `s`      - start value, if this value is 0 the values are rising with the index
/// 
/// *Note*: value is defined by (s - h * y as f64).abs(), where y is the position in the stripe
fn init(matrix: &mut ArrayViewMut2<f64>, h: f64, s: f64){
    for ((_, y), elt) in matrix.indexed_iter_mut() {
       *elt = (s - h * y as f64).abs();
    }
}

/// Initializes all four borders of both Arrays hold by the MatrixVew
/// 
/// # Arguments
/// 
/// * `matrix` - two dimensional Arrays hold by the three dimensional Matrix View
/// * `arguments` - holds the change rate h for the initialisation 
pub fn init_all(matrix: &mut ArrayViewMut3<f64>, arguments: &CalcArguments){
    init(&mut matrix.subview_mut(Axis(2),0), arguments.h, 1.0);
    init(&mut matrix.subview_mut(Axis(1),0), arguments.h, 1.0);
    init(&mut matrix.subview_mut(Axis(2),(arguments.n)), arguments.h, 0.0);
    init(&mut matrix.subview_mut(Axis(1),(arguments.n)), arguments.h, 0.0);
}

pub fn init_naive(matrix: &mut Array3<f64>, options: &Options, arguments: &CalcArguments){
    if options.inf_func == FUNC_F0 {
        let n = arguments.n;
        let h = arguments.h;
                for g in 0 .. arguments.num_matrices as usize {
                        for i in 0 .. n {
                                matrix[[g, i, 0]] = 1.0 - (h * i as f64);
                                matrix[[g, i, n]] = h * i as f64;
                                matrix[[g, 0, i]] = 1.0 - (h * i as f64);
                                matrix[[g, n, i]] = h * i as f64;
                        }

                        matrix[[g, n, 0]] = 0.0;
                        matrix[[g, 0, n]] = 0.0;
                }
        }
}



/// ONLY TESTING PURPOSE
fn init_dim1(matrix: &mut ArrayViewMut1<f64>, h: f64, s: f64){
    for (y, elt) in matrix.indexed_iter_mut() {
       *elt = (s - h * y as f64).abs();
    }
}

/// ONLY TESTING PURPOSE
#[doc(hidden)]
pub fn init_all_dim2(matrix: &mut ArrayViewMut2<f64>, arguments: &CalcArguments){
    init_dim1(&mut matrix.subview_mut(Axis(1),0), arguments.h, 1.0);
    init_dim1(&mut matrix.subview_mut(Axis(0),0), arguments.h, 1.0);
    init_dim1(&mut matrix.subview_mut(Axis(1),(arguments.n)), arguments.h, 0.0);
    init_dim1(&mut matrix.subview_mut(Axis(0),(arguments.n)), arguments.h, 0.0);
}