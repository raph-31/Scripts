use nalgebra::{DMatrix, DVector, SymmetricEigen};
pub mod turbulence;
pub use turbulence::generate_turbulent_loads;

#[derive(Debug, Clone, Copy)]
pub struct Node {
    pub x: f64,
}

impl Node {
    pub fn new(x: f64) -> Self {
        Self { x }
    }
}

pub trait Element {
    fn nodes(&self) -> (usize, usize);
    fn stiffness(&self, nodes: &[Node]) -> DMatrix<f64>;
    fn mass(&self, nodes: &[Node]) -> DMatrix<f64>;
}

#[derive(Debug, Clone, Copy)]
pub struct EulerBernoulli {
    pub n1: usize,
    pub n2: usize,
    pub e: f64,
    pub i: f64,
    pub rho: f64,
    pub a: f64,
}

impl EulerBernoulli {
    pub fn new(n1: usize, n2: usize, e: f64, i: f64, rho: f64, a: f64) -> Self {
        Self { n1, n2, e, i, rho, a }
    }
}

impl Element for EulerBernoulli {
    fn nodes(&self) -> (usize, usize) {
        (self.n1, self.n2)
    }

    fn stiffness(&self, nodes: &[Node]) -> DMatrix<f64> {
        let l = (nodes[self.n2].x - nodes[self.n1].x).abs();
        let e = self.e;
        let i = self.i;
        let k = e * i / (l.powi(3));
        let mut m = DMatrix::zeros(4, 4);
        m[(0, 0)] = 12.0 * k;
        m[(0, 1)] = 6.0 * l * k;
        m[(0, 2)] = -12.0 * k;
        m[(0, 3)] = 6.0 * l * k;

        m[(1, 0)] = 6.0 * l * k;
        m[(1, 1)] = 4.0 * l * l * k;
        m[(1, 2)] = -6.0 * l * k;
        m[(1, 3)] = 2.0 * l * l * k;

        m[(2, 0)] = -12.0 * k;
        m[(2, 1)] = -6.0 * l * k;
        m[(2, 2)] = 12.0 * k;
        m[(2, 3)] = -6.0 * l * k;

        m[(3, 0)] = 6.0 * l * k;
        m[(3, 1)] = 2.0 * l * l * k;
        m[(3, 2)] = -6.0 * l * k;
        m[(3, 3)] = 4.0 * l * l * k;
        m
    }

    fn mass(&self, nodes: &[Node]) -> DMatrix<f64> {
        let l = (nodes[self.n2].x - nodes[self.n1].x).abs();
        let m = self.rho * self.a * l / 420.0;
        let mut mm = DMatrix::zeros(4,4);
        mm[(0,0)] = 156.0 * m;
        mm[(0,1)] = 22.0 * l * m;
        mm[(0,2)] = 54.0 * m;
        mm[(0,3)] = -13.0 * l * m;

        mm[(1,0)] = 22.0 * l * m;
        mm[(1,1)] = 4.0 * l * l * m;
        mm[(1,2)] = 13.0 * l * m;
        mm[(1,3)] = -3.0 * l * l * m;

        mm[(2,0)] = 54.0 * m;
        mm[(2,1)] = 13.0 * l * m;
        mm[(2,2)] = 156.0 * m;
        mm[(2,3)] = -22.0 * l * m;

        mm[(3,0)] = -13.0 * l * m;
        mm[(3,1)] = -3.0 * l * l * m;
        mm[(3,2)] = -22.0 * l * m;
        mm[(3,3)] = 4.0 * l * l * m;
        mm
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Timoshenko {
    pub n1: usize,
    pub n2: usize,
    pub e: f64,
    pub i: f64,
    pub g: f64,
    pub a: f64,
    pub kappa: f64,
    pub rho: f64,
}

impl Timoshenko {
    pub fn new(n1: usize, n2: usize, e: f64, i: f64, g: f64, a: f64, kappa: f64, rho: f64) -> Self {
        Self { n1, n2, e, i, g, a, kappa, rho }
    }
}

impl Element for Timoshenko {
    fn nodes(&self) -> (usize, usize) {
        (self.n1, self.n2)
    }

    fn stiffness(&self, nodes: &[Node]) -> DMatrix<f64> {
        let l = (nodes[self.n2].x - nodes[self.n1].x).abs();
        let e = self.e;
        let i = self.i;
        let g = self.g;
        let a = self.a;
        let kappa = self.kappa;
        let phi = 12.0 * e * i / (kappa * a * g * l * l);
        let k = e * i / (l.powi(3) * (1.0 + phi));

        let mut m = DMatrix::zeros(4, 4);
        m[(0, 0)] = 12.0 * k;
        m[(0, 1)] = 6.0 * l * k;
        m[(0, 2)] = -12.0 * k;
        m[(0, 3)] = 6.0 * l * k;

        m[(1, 0)] = 6.0 * l * k;
        m[(1, 1)] = (4.0 + phi) * l * l * k;
        m[(1, 2)] = -6.0 * l * k;
        m[(1, 3)] = (2.0 - phi) * l * l * k;

        m[(2, 0)] = -12.0 * k;
        m[(2, 1)] = -6.0 * l * k;
        m[(2, 2)] = 12.0 * k;
        m[(2, 3)] = -6.0 * l * k;

        m[(3, 0)] = 6.0 * l * k;
        m[(3, 1)] = (2.0 - phi) * l * l * k;
        m[(3, 2)] = -6.0 * l * k;
        m[(3, 3)] = (4.0 + phi) * l * l * k;
        m
    }

    fn mass(&self, nodes: &[Node]) -> DMatrix<f64> {
        let l = (nodes[self.n2].x - nodes[self.n1].x).abs();
        let m = self.rho * self.a * l / 420.0;
        let mut mm = DMatrix::zeros(4,4);
        mm[(0,0)] = 156.0 * m;
        mm[(0,1)] = 22.0 * l * m;
        mm[(0,2)] = 54.0 * m;
        mm[(0,3)] = -13.0 * l * m;

        mm[(1,0)] = 22.0 * l * m;
        mm[(1,1)] = 4.0 * l * l * m;
        mm[(1,2)] = 13.0 * l * m;
        mm[(1,3)] = -3.0 * l * l * m;

        mm[(2,0)] = 54.0 * m;
        mm[(2,1)] = 13.0 * l * m;
        mm[(2,2)] = 156.0 * m;
        mm[(2,3)] = -22.0 * l * m;

        mm[(3,0)] = -13.0 * l * m;
        mm[(3,1)] = -3.0 * l * l * m;
        mm[(3,2)] = -22.0 * l * m;
        mm[(3,3)] = 4.0 * l * l * m;
        mm
    }
}

pub enum BeamElement {
    Euler(EulerBernoulli),
    Timo(Timoshenko),
}

impl Element for BeamElement {
    fn nodes(&self) -> (usize, usize) {
        match self {
            BeamElement::Euler(e) => e.nodes(),
            BeamElement::Timo(t) => t.nodes(),
        }
    }

    fn stiffness(&self, nodes: &[Node]) -> DMatrix<f64> {
        match self {
            BeamElement::Euler(e) => e.stiffness(nodes),
            BeamElement::Timo(t) => t.stiffness(nodes),
        }
    }

    fn mass(&self, nodes: &[Node]) -> DMatrix<f64> {
        match self {
            BeamElement::Euler(e) => e.mass(nodes),
            BeamElement::Timo(t) => t.mass(nodes),
        }
    }
}

pub struct Beam {
    pub nodes: Vec<Node>,
    pub elements: Vec<BeamElement>,
}

impl Beam {
    pub fn new() -> Self {
        Self { nodes: Vec::new(), elements: Vec::new() }
    }

    pub fn add_node(&mut self, x: f64) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(Node::new(x));
        idx
    }

    pub fn add_element(&mut self, elem: BeamElement) {
        self.elements.push(elem);
    }

    pub fn assemble(&self) -> DMatrix<f64> {
        let ndof = self.nodes.len() * 2;
        let mut k_global = DMatrix::<f64>::zeros(ndof, ndof);
        for elem in &self.elements {
            let (i, j) = elem.nodes();
            let k_local = elem.stiffness(&self.nodes);
            let map = [2*i, 2*i+1, 2*j, 2*j+1];
            for a in 0..4 {
                for b in 0..4 {
                    k_global[(map[a], map[b])] += k_local[(a, b)];
                }
            }
        }
        k_global
    }

    pub fn assemble_mass(&self) -> DMatrix<f64> {
        let ndof = self.nodes.len() * 2;
        let mut m_global = DMatrix::<f64>::zeros(ndof, ndof);
        for elem in &self.elements {
            let (i, j) = elem.nodes();
            let m_local = elem.mass(&self.nodes);
            let map = [2*i, 2*i+1, 2*j, 2*j+1];
            for a in 0..4 {
                for b in 0..4 {
                    m_global[(map[a], map[b])] += m_local[(a, b)];
                }
            }
        }
        m_global
    }

    pub fn natural_frequencies(&self, modes: usize) -> Vec<f64> {
        let k = self.assemble();
        let m = self.assemble_mass();
        let k_red = k.slice_range(2.., 2..).into_owned();
        let m_red = m.slice_range(2.., 2..).into_owned();
        let m_chol = m_red.cholesky().unwrap();
        let m_inv = m_chol.inverse();
        let a = &m_inv.transpose() * k_red * &m_inv;
        let eig = SymmetricEigen::new(a);
        let mut freqs: Vec<f64> = eig
            .eigenvalues
            .iter()
            .map(|&l| (l.max(0.0)).sqrt() / (2.0 * std::f64::consts::PI))
            .collect();
        freqs.sort_by(|a, b| a.partial_cmp(b).unwrap());
        freqs.truncate(modes);
        freqs
    }

    pub fn solve(&self, loads: &[f64]) -> DVector<f64> {
        let k_global = self.assemble();
        let f = DVector::from_column_slice(loads);
        // For simplicity assume first node fixed (w=0, theta=0)
        let mut k_reduced = k_global.slice_range(2.., 2..).into_owned();
        let f_reduced = f.rows(2, f.len() - 2).into_owned();
        let disp_reduced = k_reduced.lu().solve(&f_reduced).unwrap();
        let mut disp = DVector::zeros(f.len());
        disp.rows_mut(2, disp_reduced.len()).copy_from(&disp_reduced);
        disp
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euler_beam() {
        let mut beam = Beam::new();
        beam.add_node(0.0);
        beam.add_node(1.0);
        beam.add_element(BeamElement::Euler(EulerBernoulli::new(0,1,210e9,8.33e-6,7850.0,1e-2)));
        let f = vec![0.0,0.0,1000.0,0.0];
        let d = beam.solve(&f);
        assert!(d[2].abs() > 0.0);
    }

    #[test]
    fn timoshenko_beam() {
        let mut beam = Beam::new();
        beam.add_node(0.0);
        beam.add_node(1.0);
        beam.add_element(BeamElement::Timo(Timoshenko::new(0,1,210e9,8.33e-6,80e9,1e-2,5.0/6.0,7850.0)));
        let f = vec![0.0,0.0,1000.0,0.0];
        let d = beam.solve(&f);
        assert!(d[2].abs() > 0.0);
    }

    #[test]
    fn modal_comparison() {
        let mut beam_e = Beam::new();
        beam_e.add_node(0.0);
        beam_e.add_node(1.0);
        beam_e.add_element(BeamElement::Euler(EulerBernoulli::new(0,1,210e9,8.33e-6,7850.0,1e-2)));
        let fe = beam_e.natural_frequencies(1)[0];

        let mut beam_t = Beam::new();
        beam_t.add_node(0.0);
        beam_t.add_node(1.0);
        beam_t.add_element(BeamElement::Timo(Timoshenko::new(0,1,210e9,8.33e-6,80e9,1e-2,5.0/6.0,7850.0)));
        let ft = beam_t.natural_frequencies(1)[0];

        assert!((fe - ft).abs() / fe < 0.2);
    }
}
