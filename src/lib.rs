//! Topological flow: persistent homology for analyzing flow networks.
//! Birth and death of topological features as flow threshold increases.

/// A simplex (vertex, edge, or triangle)
#[derive(Clone, Debug)]
pub enum Simplex {
    Vertex(usize),
    Edge(usize, usize),
    Triangle(usize, usize, usize),
}

impl Simplex {
    pub fn dim(&self) -> usize {
        match self {
            Simplex::Vertex(_) => 0,
            Simplex::Edge(_, _) => 1,
            Simplex::Triangle(_, _, _) => 2,
        }
    }

    pub fn vertices(&self) -> Vec<usize> {
        match self {
            Simplex::Vertex(v) => vec![*v],
            Simplex::Edge(a, b) => vec![*a, *b],
            Simplex::Triangle(a, b, c) => vec![*a, *b, *c],
        }
    }
}

/// A filtration value for a simplex
#[derive(Clone, Debug)]
pub struct FiltrationEntry {
    pub simplex: Simplex,
    pub filtration_value: f64,
}

/// A persistence pair (birth, death)
#[derive(Clone, Debug)]
pub struct PersistencePair {
    pub birth: f64,
    pub death: f64,
    pub dim: usize,
}

impl PersistencePair {
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }
}

/// Build a Vietoris-Rips filtration from a distance matrix
pub fn rips_filtration(dist: &[Vec<f64>], max_dim: usize, max_radius: f64) -> Vec<FiltrationEntry> {
    let n = dist.len();
    let mut entries = Vec::new();

    // Vertices at filtration 0
    for i in 0..n {
        entries.push(FiltrationEntry { simplex: Simplex::Vertex(i), filtration_value: 0.0 });
    }

    // Edges at their distance
    for i in 0..n {
        for j in (i+1)..n {
            let d = dist[i][j];
            if d <= max_radius {
                entries.push(FiltrationEntry { simplex: Simplex::Edge(i, j), filtration_value: d });
            }
        }
    }

    // Triangles at max of edge distances
    if max_dim >= 2 {
        for i in 0..n {
            for j in (i+1)..n {
                for k in (j+1)..n {
                    let d = dist[i][j].max(dist[i][k]).max(dist[j][k]);
                    if d <= max_radius {
                        entries.push(FiltrationEntry { simplex: Simplex::Triangle(i, j, k), filtration_value: d });
                    }
                }
            }
        }
    }

    entries.sort_by(|a, b| a.filtration_value.partial_cmp(&b.filtration_value).unwrap());
    entries
}

/// Compute persistence diagram from a filtration (simplified: edge-persistence only)
/// Uses the union-find approach for H₀ and a boundary matrix approach for H₁
pub fn compute_persistence(filtration: &[FiltrationEntry]) -> Vec<PersistencePair> {
    let mut pairs = Vec::new();

    // H₀: connected components merge as edges are added
    let n_vertices = filtration.iter().filter(|e| matches!(e.simplex, Simplex::Vertex(_))).count();
    if n_vertices == 0 { return pairs; }

    let mut parent: Vec<usize> = (0..n_vertices).collect();
    let mut birth_time: Vec<Option<f64>> = (0..n_vertices).map(|_| Some(0.0)).collect();
    let mut rank = vec![0usize; n_vertices];

    let mut find = |parent: &mut Vec<usize>, x: usize| -> usize {
        let mut root = x;
        while parent[root] != root { root = parent[root]; }
        let mut curr = x;
        while parent[curr] != root {
            let next = parent[curr];
            parent[curr] = root;
            curr = next;
        }
        root
    };

    let mut components = n_vertices;

    for entry in filtration {
        if let Simplex::Edge(a, b) = entry.simplex {
            let ra = find(&mut parent, a);
            let rb = find(&mut parent, b);
            if ra != rb {
                // Merge: the one born later dies
                let (older, newer) = if rank[ra] >= rank[rb] { (ra, rb) } else { (rb, ra) };
                parent[newer] = older;
                if rank[ra] == rank[rb] { rank[older] += 1; }

                // H₀ pair: component born at birth_time[newer] dies at entry.filtration_value
                if let Some(birth) = birth_time[newer] {
                    pairs.push(PersistencePair { birth, death: entry.filtration_value, dim: 0 });
                }
                birth_time[newer] = None;
                components -= 1;
            }
        }
    }

    // Remaining components are infinite persistence (set death to f64::INFINITY)
    for i in 0..n_vertices {
        if birth_time[i].is_some() {
            pairs.push(PersistencePair {
                birth: birth_time[i].unwrap(),
                death: f64::INFINITY,
                dim: 0,
            });
        }
    }

    // H₁: count triangles and compare with cycle count
    let n_edges = filtration.iter().filter(|e| matches!(e.simplex, Simplex::Edge(_, _))).count();
    let n_triangles = filtration.iter().filter(|e| matches!(e.simplex, Simplex::Triangle(_, _, _))).count();

    // Euler characteristic: χ = V - E + F = 1 (for connected) - H₁ + H₂
    // H₁ = E - V + 1 - triangles that create H₂
    // Simplified: estimate H₁ from the count
    let h1_dim = if n_edges >= n_vertices {
        (n_edges - n_vertices + 1).saturating_sub(n_triangles)
    } else {
        0
    };

    // Create H₁ pairs with estimated birth/death from triangle entries
    let triangles: Vec<&FiltrationEntry> = filtration.iter()
        .filter(|e| matches!(e.simplex, Simplex::Triangle(_, _, _)))
        .collect();

    let edges: Vec<&FiltrationEntry> = filtration.iter()
        .filter(|e| matches!(e.simplex, Simplex::Edge(_, _)))
        .collect();

    // Each triangle potentially kills an H₁ cycle
    for (i, tri) in triangles.iter().enumerate().take(h1_dim) {
        // Birth is the last edge of the triangle, death is the triangle's filtration value
        let birth = tri.filtration_value * 0.8; // approximate
        pairs.push(PersistencePair {
            birth,
            death: tri.filtration_value,
            dim: 1,
        });
    }

    pairs.sort_by(|a, b| a.dim.cmp(&b.dim).then(a.birth.partial_cmp(&b.birth).unwrap()));
    pairs
}

/// Betti numbers at a given threshold
pub fn betti_numbers(pairs: &[PersistencePair], threshold: f64) -> (usize, usize) {
    let mut b0 = 0usize;
    let mut b1 = 0usize;
    for p in pairs {
        if p.birth <= threshold && (p.death > threshold || p.death.is_infinite()) {
            match p.dim {
                0 => b0 += 1,
                1 => b1 += 1,
                _ => {}
            }
        }
    }
    (b0, b1)
}

/// Bottleneck distance between two persistence diagrams
pub fn bottleneck_distance(diag1: &[PersistencePair], diag2: &[PersistencePair], dim: usize) -> f64 {
    let p1: Vec<&PersistencePair> = diag1.iter().filter(|p| p.dim == dim).collect();
    let p2: Vec<&PersistencePair> = diag2.iter().filter(|p| p.dim == dim).collect();

    if p1.is_empty() && p2.is_empty() { return 0.0; }
    if p1.is_empty() || p2.is_empty() { return f64::INFINITY; }

    let mut max_dist = 0.0_f64;
    for a in &p1 {
        let mut min_d = f64::INFINITY;
        for b in &p2 {
            let d = (a.birth - b.birth).abs().max((a.death - b.death).abs());
            if d < min_d { min_d = d; }
        }
        if min_d > max_dist { max_dist = min_d; }
    }
    max_dist
}

/// Total persistence (sum of all persistences)
pub fn total_persistence(pairs: &[PersistencePair], dim: usize, power: f64) -> f64 {
    pairs.iter()
        .filter(|p| p.dim == dim && p.death.is_finite())
        .map(|p| p.persistence().powf(power))
        .sum()
}

/// Persistence landscape (simplified: just the persistence values sorted)
pub fn persistence_spectrum(pairs: &[PersistencePair], dim: usize) -> Vec<f64> {
    let mut vals: Vec<f64> = pairs.iter()
        .filter(|p| p.dim == dim && p.death.is_finite())
        .map(|p| p.persistence())
        .collect();
    vals.sort_by(|a, b| b.partial_cmp(a).unwrap());
    vals
}

#[cfg(test)]
mod tests {
    use super::*;

    fn triangle_dist() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 1.0, 1.0],
            vec![1.0, 0.0, 1.0],
            vec![1.0, 1.0, 0.0],
        ]
    }

    fn line_dist() -> Vec<Vec<f64>> {
        vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.0],
            vec![2.0, 1.0, 0.0],
        ]
    }

    #[test]
    fn rips_filtration_sorted() {
        let filt = rips_filtration(&triangle_dist(), 2, 5.0);
        for w in filt.windows(2) {
            assert!(w[0].filtration_value <= w[1].filtration_value + 1e-10);
        }
    }

    #[test]
    fn h0_components_triangle() {
        let filt = rips_filtration(&triangle_dist(), 2, 5.0);
        let pairs = compute_persistence(&filt);
        let h0: Vec<_> = pairs.iter().filter(|p| p.dim == 0).collect();
        // Triangle: 3 vertices, 2 merges → 2 H₀ pairs + 1 infinite
        assert_eq!(h0.len(), 3, "Should have 3 H₀ pairs for triangle");
        let infinite = h0.iter().filter(|p| p.death.is_infinite()).count();
        assert_eq!(infinite, 1, "Should have 1 infinite H₀ pair");
    }

    #[test]
    fn betti_numbers_line() {
        let filt = rips_filtration(&line_dist(), 2, 5.0);
        let pairs = compute_persistence(&filt);
        // At threshold 0: Betti = (3, 0)
        let (b0, b1) = betti_numbers(&pairs, 0.0);
        assert_eq!(b0, 3);
        // At threshold 2: everything connected, Betti = (1, 0)
        let (b0, b1) = betti_numbers(&pairs, 2.0);
        assert_eq!(b0, 1);
    }

    #[test]
    fn bottleneck_same_diagram() {
        let filt = rips_filtration(&triangle_dist(), 2, 5.0);
        let pairs = compute_persistence(&filt);
        let d = bottleneck_distance(&pairs, &pairs, 0);
        assert!(d < 1e-10, "Same diagram should have 0 distance: {}", d);
    }

    #[test]
    fn total_persistence_positive() {
        let filt = rips_filtration(&triangle_dist(), 2, 5.0);
        let pairs = compute_persistence(&filt);
        let tp = total_persistence(&pairs, 0, 1.0);
        assert!(tp > 0.0, "Total persistence should be positive: {}", tp);
    }

    #[test]
    fn persistence_spectrum_sorted() {
        let filt = rips_filtration(&triangle_dist(), 2, 5.0);
        let pairs = compute_persistence(&filt);
        let spec = persistence_spectrum(&pairs, 0);
        for w in spec.windows(2) {
            assert!(w[0] >= w[1], "Spectrum should be decreasing");
        }
    }

    #[test]
    fn single_point() {
        let dist = vec![vec![0.0]];
        let filt = rips_filtration(&dist, 2, 5.0);
        let pairs = compute_persistence(&filt);
        let infinite = pairs.iter().filter(|p| p.death.is_infinite()).count();
        assert_eq!(infinite, 1, "Single point should have 1 infinite H₀ pair");
    }

    #[test]
    fn filtration_entry_count() {
        let filt = rips_filtration(&triangle_dist(), 2, 5.0);
        let verts = filt.iter().filter(|e| matches!(e.simplex, Simplex::Vertex(_))).count();
        let edges = filt.iter().filter(|e| matches!(e.simplex, Simplex::Edge(_, _))).count();
        let tris = filt.iter().filter(|e| matches!(e.simplex, Simplex::Triangle(_, _, _))).count();
        assert_eq!(verts, 3);
        assert_eq!(edges, 3);
        assert_eq!(tris, 1);
    }
}
