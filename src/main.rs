use std::time::Instant;

use alco::Octavian;
use std::collections::HashSet;



fn main() {
    println!("Hello world");  
    // Start the timer
    let start = Instant::now();

    let b = Octavian::<isize>::BASIS;

    for x in b {
        println!("{:?}", x*x);
        println!("{:?}", x.conj());
    }

    // // Your main logic here
    // let mut units: HashSet<Octavian<i8>> = HashSet::new();
    // for u in Octavian::<i8>::OCTAVIAN_UNITS_COEFFICIENTS {
    //     let x = Octavian::new(u);
    //     units.insert(x);
    // }
    // assert_eq!(240, units.len());
    // let mut result = HashSet::<Octavian<i8>>::new();
    // for u in &units {
    //     for v in &units {
    //         result.insert(u.clone() * v.clone());
    //     }
    // }
    // assert_eq!(240, result.len());

    // for x in OCTAVIAN_UNITS_COEFFICIENTS {
    //     let y = Octavian::new(x);
    //     println!("{:?}", y.norm()); 
    // }

    // let y = Octavian::<i64>::one();
    // let y = Octavian::<i32>::new([ 1,2,3,4,5,6,7,8 ]);
    // let y = y*y;
    // println!("{:?}", y);
    // for row in y.left_adjoint_matrix() {
    //     println!("{:?}", row);
    // }

    // println!("{:?}", Octavian::<i32>::one().scale(y.trace())); 

    // Stop the timer and print the elapsed time
    let duration = start.elapsed();
    println!("Execution time: {:?}", duration);
}

