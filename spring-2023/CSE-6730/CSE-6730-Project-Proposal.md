# Integrating Audio Generation into 2D Artificial Life Models

There are several visually stunning, and sometimes convincingly lifelike, examples of artificial life. [1][2][3]

As an audio rendering approach, I propose rendering measurements from vertices of convex hulls of each (visually) rendered ellipse.

Measurable `Node` values:

- $Position \in \mathbb{R}^{N}, N \in (2, 3)$
- Measurable values from any of $N$-nearest neighbor vertices
- Incoming attractive force from any (or all, or any linear combination of) other vertices
- Outgoing repulsive force to any (or all, or any linear combination of) other vertices

Render to audio in real-time, with configurable rendering properties:

- $SR_m$: **Measurement sample rate**
  - How often to sample the measured values of the (ellipsoid convex hull) vertices of the monitored nodes, $\mathbf{Nodes}_\text{Monitored}$.
- $\mathbf{Nodes}_\text{Monitored}$: **Monitored nodes**
  - All properties on each instance hold most recently measured values
    - Sampled at configurable rate $SR_m \in \mathbb{I}$ (measurement sample frequency in $1/s$).
- $f(\mathbf{nodes}), \mathbf{nodes} \in \mathbb{Node}^{N}, N \in \mathbb{I}$:

## Approach

- Start from an existing artificial life implementation, like Tom Mohr's `particle-life-app`[1].
- Note that it's entirely composed of tiny filled ellipsis of different color (and look to be all circles even?)
- Proposed model additions for audio generation:
  - (See [C++ Definitions[#c++-definitions] at the bottom of this document.)
  - Render 2D ellipsoids to audio
    - Keep track of the $N$-sided 2D convex hull surrounding each ellipse.
      - As $N$ increases, the convex hull approaches the ellipse (more granular).
      - Add a slider to control the number of sides.
    - User can "tap" the current attractional force between any two vertices of the convex hulls surrounding any ellipse.
    - Let $\mathbb{N}$ be the set of all convex hull vertices.
    - User can attach a sensor to any convex hull vertex, and measure its:
      - Nearest $n$ nodes (along with any of their measured properties)
      - `position`
        - `position` derivative (available as `const public &velocity`)
      - Total incoming force (how much is this node being attracted to any set of monitored nodes)
        - Or any combination of point's measured, e.g. total incoming force (normalized to total incoming force across all nodes):
          - $\dfrac{\sum\limits_{\mathbf{n} \in \mathbf{Neighbors}}\mathbf{F_{In}}\left(\mathbf{n}\right)}{\sum\limits_{\mathbf{n} \in \mathbf{Nodes}}\mathbf{F_{In}}\left(\mathbf{n}\right)}$
  - Extend 2D ellipse into a 3D ellipsoid.
    - Add a toggles to existing ImGui app to:
      - Switch from 2D ellipse visual rendering to 3D
      - Switch audio rendering to be based on 3D ellipsoid rather than 2D ellipse
    - Add a toggle to existing app to switch from current 2D ellipse rendering to 3D. Render the 3D ellipsoids to audio

## Citations

[1] T. Mohr, “🦠 Particle Life App.” Feb. 04, 2023. Accessed: Feb. 03, 2023. [Online]. Available: https://github.com/tom-mohr/particle-life-app
[2] TODO
[3] TODO

## C++ definitions

```cpp
struct Value:
    const double &value = Value = Value; // Current value
    const vector<double> = Deltas = Deltas; // $v_t - v_{t-1}`$

    // Stored from most-to-least recent.
    // Length configurable.
    // Defaults to the single-element history `{most_recent_value}`.
    const vector<double> &ValueHistory = LatestValues;

    const Set(const ValueType new_value) {
        ValueHistory.emplace_back(new_value);
        Deltas.emplace_back(new_value - Value);
        Value = new_value;
    }

    static const long DefaultMeasureSampleRate = 48_000; // Sampling frequency in 1/s

private:
    double Value;
    vector<double> Deltas;
    vector<double> ValueHistory;
```

```cpp
struct AudioRenderer: {
    const vector<Node> &MonitoredNodes;
    // How often to sample the vertices of the convex hulls of `MonitoredNodes`.
    const long MeasureSampleRate = &MeasureSampleRate;

private:
    long MeasureSampleRate = DefaultMeasureSampleRate;
}
```
