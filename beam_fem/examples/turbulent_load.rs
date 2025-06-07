use beam_fem::turbulence::generate_turbulent_loads;

fn main() {
    let loads = generate_turbulent_loads(1.0, 2.0, 5.0, 0.05, 3, 20);
    for (i, step) in loads.iter().enumerate().take(3) {
        println!("step {i}: {step:?}");
    }
}
