use beam_fem::{BeamElement, Beam, EulerBernoulli, Timoshenko};

fn main() {
    let mut beam_e = Beam::new();
    beam_e.add_node(0.0);
    beam_e.add_node(1.0);
    beam_e.add_element(BeamElement::Euler(EulerBernoulli::new(0, 1, 210e9, 8.33e-6, 7850.0, 1e-2)));
    let freqs_e = beam_e.natural_frequencies(2);
    println!("Euler-Bernoulli first two frequencies: {:?}", freqs_e);

    let mut beam_t = Beam::new();
    beam_t.add_node(0.0);
    beam_t.add_node(1.0);
    beam_t.add_element(BeamElement::Timo(Timoshenko::new(0, 1, 210e9, 8.33e-6, 80e9, 1e-2, 5.0/6.0, 7850.0)));
    let freqs_t = beam_t.natural_frequencies(2);
    println!("Timoshenko first two frequencies: {:?}", freqs_t);
}
