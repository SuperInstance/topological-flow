# topological-flow

**Persistent homology for flow networks — track the birth and death of topological features as flow increases.**

Pure Rust, zero dependencies. Builds Vietoris-Rips filtrations from distance matrices and computes persistence diagrams that reveal the topological structure of flow networks at every scale.

## What This Gives You

- **Vietoris-Rips filtration** — simplices appear as the distance threshold increases
- **Persistence diagrams** — birth/death pairs for 0-dimensional (components) and 1-dimensional (loops) features
- **Betti numbers** — count of topological features at each threshold
- **Zero dependencies** — all from scratch

## The Core Idea

A flow network has structure at multiple scales. At low flow, the network is fragmented — many isolated components. As flow increases, components merge, loops form, and the topology changes. Persistent homology tracks *which* features survive across scales and *how long* they persist.

A feature that's born early and dies late is a robust topological property. A feature that appears and immediately vanishes is noise.

## Quick Start

```rust
use topological_flow::{rips_filtration, compute_persistence, PersistencePair};

// Distance matrix for 4 points
let dist = vec![
    vec![0.0, 1.0, 2.0, 3.0],
    vec![1.0, 0.0, 1.0, 2.0],
    vec![2.0, 1.0, 0.0, 1.0],
    vec![3.0, 2.0, 1.0, 0.0],
];

// Build filtration up to radius 3.0, allowing triangles
let filtration = rips_filtration(&dist, 2, 3.0);

// Compute persistence pairs
let pairs = compute_persistence(&filtration);
for pair in &pairs {
    println!("H{}: born at {:.2}, dies at {:.2} (persistence = {:.2})",
        pair.dim, pair.birth, pair.death, pair.persistence());
}
```

## API Reference

| Type / Function | Description |
|----------------|-------------|
| `Simplex` | Vertex, Edge, or Triangle |
| `FiltrationEntry` | A simplex with its filtration value |
| `PersistencePair` | Birth, death, and dimension of a topological feature |
| `rips_filtration(dist, max_dim, max_radius)` | Build Vietoris-Rips filtration |
| `compute_persistence(filtration)` | Compute persistence diagram |

## How It Fits

Part of the SuperInstance topological ecosystem:

- **[topology-lab](https://github.com/SuperInstance/topology-lab)** — Interactive visualization of persistent homology
- **[wasserstein-narrative](https://github.com/SuperInstance/wasserstein-narrative)** — Persistence diagrams for story analysis
- **topological-flow** — Persistence for flow networks (this repo)

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
topological-flow = { git = "https://github.com/SuperInstance/topological-flow" }
```

## License

MIT

Part of the [SuperInstance](https://github.com/SuperInstance) ecosystem.
