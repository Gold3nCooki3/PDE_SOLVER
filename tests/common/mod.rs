#![allow(unused)]
extern crate pde_solver;
extern crate ndarray;

use self::pde_solver::basis::askparams;
use self::pde_solver::basis::displaymatrix;
use self::pde_solver::calculate::{serial,parallel};
use self::pde_solver::basis::*;
use self::ndarray::arr3;