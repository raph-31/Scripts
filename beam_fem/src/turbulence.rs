//! Simple utilities to generate turbulent line forces.
//!
//! The generated field is a random process with optional spatial
//! smoothing (controlled by the correlation length) and a basic
//! temporal filter to limit the highest frequency content.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};

/// Generate a time history of turbulent line forces.
///
/// * `amplitude` - standard deviation of the force per node
/// * `correlation_length` - spatial correlation length in number of nodes
/// * `freq_cutoff` - upper frequency cutoff in Hz for a first-order low pass
/// * `dt` - time step in seconds
/// * `nodes` - number of nodes along the beam
/// * `steps` - number of time steps
pub fn generate_turbulent_loads(
    amplitude: f64,
    correlation_length: f64,
    freq_cutoff: f64,
    dt: f64,
    nodes: usize,
    steps: usize,
) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(0);
    let normal = Normal::new(0.0, amplitude).unwrap();
    let mut data = vec![vec![0.0; nodes]; steps];

    // spatial smoothing window size (in nodes)
    let corr = correlation_length.max(0.0).round() as usize;

    // coefficient for simple exponential low-pass filter
    let alpha = (-2.0 * std::f64::consts::PI * freq_cutoff * dt).exp();

    // generate uncorrelated noise first
    for t in 0..steps {
        for n in 0..nodes {
            data[t][n] = normal.sample(&mut rng);
        }
    }

    // spatial smoothing using moving average
    if corr > 1 {
        for t in 0..steps {
            let mut smoothed = vec![0.0; nodes];
            for i in 0..nodes {
                let start = i.saturating_sub(corr);
                let end = (i + corr + 1).min(nodes);
                smoothed[i] = data[t][start..end].iter().sum::<f64>() / (end - start) as f64;
            }
            data[t] = smoothed;
        }
    }

    // temporal low-pass filter
    if freq_cutoff > 0.0 {
        for n in 0..nodes {
            let mut state = 0.0;
            for t in 0..steps {
                state = alpha * state + (1.0 - alpha) * data[t][n];
                data[t][n] = state;
            }
        }
    }

    data
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape() {
        let out = generate_turbulent_loads(1.0, 2.0, 5.0, 0.1, 4, 10);
        assert_eq!(out.len(), 10);
        assert_eq!(out[0].len(), 4);
    }
}

