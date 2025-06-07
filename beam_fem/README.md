# beam_fem

A lightweight Rust library for 1D beam finite element analysis, inspired by the Python library `beef`.
It supports Euler-Bernoulli and Timoshenko beam elements.

```
use beam_fem::{BeamElement, Beam, EulerBernoulli, Timoshenko};

fn main() {
    let mut beam = Beam::new();
    beam.add_node(0.0);
    beam.add_node(1.0);

    let e = 210e9;        // Young's modulus
    let i = 8.33e-6;      // Second moment of area
    let a = 1e-2;         // Area
    let g = 80e9;         // Shear modulus
    let rho = 7850.0;     // Density

    beam.add_element(BeamElement::Euler(EulerBernoulli::new(0, 1, e, i, rho, a)));
    let loads = vec![0.0, 0.0, 1000.0, 0.0]; // nodal load at node 2
    let disp = beam.solve(&loads);
    println!("{:?}", disp);

    let freqs = beam.natural_frequencies(1);
    println!("first mode: {} Hz", freqs[0]);
}
```

See `examples` directory for more usage details.

### Modal analysis example

To compute the first two natural frequencies for Euler-Bernoulli and
Timoshenko beams run:

```bash
cargo run --example modal
```
