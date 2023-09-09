# Notes for 8803

**Current plan:**
- I'm _committing_ to adding damping estimation to mesh2audio, and I'm punting on the low-rank updates.
- I'm _committing_ to near-real-time re-estimation of mesh2audio model via fast eigenvalue re-estimation by restricting vertex position changes to low-rank matrix updates.
    - The problem I'm solving is that generating many synthetic training samples is going to be prohibitively slowe if we need to do fresh eigenvalue estimation each time.
    - At least at first - only support random walks along (low-rank position delta, material properties) space.
    - Later we can support finding the low-rank delta resulting in positions nearest to a target mesh.

## Project ideas

$\Phi$ = The generative distribution representing my dream synthetic, realtime-generated ✨mesh+audio+more✨ dataset.

A single sample $\phi \sim \Phi$ has the following attributes:
- Mesh
- Material params
- Excitation vertices
- Excitation "hammer" audio model params for each vertex
- Audio
    - Audio is generated via a source-filter modal synthesis modal.
- Video
    - Video is generated via physically-based rendering (PBR) of the mesh and material properties.

A stream of $[\phi_0, \cdots, \phi_t] \sim \Phi^{T}$ can be generated (in real time?!) using any of the following methods:
- Known (mesh+excitation vertices) stream.
    - Implement an automatic generator using the following procedure:
      - Given a set of known meshes (sampled from any big mesh library), 
      - Remain "object-like".
- Video
    - Estimate meshes and excitation vertices from a video stream.
- Video+Audio (**This is a potential followup project.**)
    - Improving estimation of mesh/material/interaction vertices streams by incorporating audio data.


I want to extend my mesh2audio work.
I have a good path forward for my broader research project.
It would be nice to keep this a linear path, and do this stuff before my MLG project, as prerequisites.

I'm not even positive there's a very obvious GNN tie-in after getting the higher-quality (via damping matrix estimation) and faster mesh2audio.
But I have a strong intuition that being able to generate, in near-real-time, good quality synthetic $\phi \sim \Phi$ samples that are with fully-known "true" parameters, will unlock many potential GNN applications.

Why not decouple the MLG project from all mesh2audio work?
    - I'll probably never do it otherwise 🤷 😬
    - I really should be spending lots of time on my research.
      Grad school isn't just about taking classes and getting good grades!
      It would be very nice to have something strong to show for my two years in grad school.

Main problem: Need to re-solve the eigenvalue problem every time the graph changes, making it slow to generate audio samples from meshes for an accuracy reading.
  - Is there a way to more cheaply estimate an eigenvalue update based on a known matrix/eigens and a matrix update?
    - Guessing it would have to be a low-rank matrix update.
      Well, for generating synthetic datasets of $\Phi$ we could enforce that all mesh position delta matrices are low-rank.
      This would allow for fast regeneration of the (Faust) audio model, potentially making it feasible to generate very large synthetic data corpuses for self-supervised learning (or to augment existing video+audio datasets with physically plausible generated samples of (audio, mesh, material props, excitation positions,modal synthesis props (including "hammer" params, estimated modes, damping, ...)
    - Blocker honestly might be getting mesh2audio to estimate damping in the finite element analysis.


Yup! That's it: Next big step is to including damping estimation in mesh2audio (and mesh2faust).
  - I think I should treat this as a blocker for my MLG project.
    It's something that really needs to be done for mesh2audio, and there's not really another mesh2audio out there.
    Think how cool it would be to have the synthetic dataset (see "Main problem" above for specifics) with _completely_ known parameters.
    For use in self/semi-supervised learning, or audio+video dataset augmentation.
    It opens up several projects listed below, e.g. improving 3d mesh generation from video by including audio.
    Combine with fast synthetic data generation via fast eigenvalue-reestimation with low-rank mesh position deltas.

## Next steps

These next steps apply more broadly than just MLG.
But I want to treat these as blockers for MLG work, and front-load them.

- Merge mesh2faust work.
    - Post on GH issue about not being able to reach Romain.
    - Mention my interest in extending the work to include damping matrix estimation.
        - (Get the wording right on damping estimation - look up in FEM literature.)
- Duplicate mesh2audio repo to github/khiner/mesh2audio.
    - Private for now. Link to the current (class project) repo and don't touch it anymore.
- Add damping estimation to mesh2audio.
- Near-real-time re-estimation of mesh2audio model via fast eigenvalue re-estimation after changing mesh positions.
    - Investigate fast eigenvalue re-estimation after mesh position updates.
    - Restrict vertex position changes to low-rank matrix deltas.
        - Using low-rank eigenvalue update techniques
    - Later (stretch goal):
        - Given an arbitrary mesh target, find the closest low-rank position update.
            - We already have the eigenvalues and eigenvectors for the current mesh positions.
            - Repeat until convergence:
                - Find the _low-rank_ position delta that produces the mesh _most similar to_ the target mesh. (Minimize $d(mesh_\text{current}, mesh_\{target}$)).
        - In mesh2audio UI: Show a dim overlay of the computed closest low-rank updated positions used for the current "fast" audio model update based on the low-rank eigenvalue update.
        - Maybe: Kick off a full eigenvalue re-estimation job in the background. (Low-rank update is just a fast audio re-esimation.)
            - Only useful in the app context - not necessary for fast synthetic dataset generation, since we can generate the (mesh, material props) samples however we want if there's no target mesh, e.g. a random walk over low-rank position updates.

Potential projects
- Maybe: Learn a mesh2audio model by _directly_ mapping the graph to the audio samples.
  - This is me trying to be as GNN-forward as possible.
- Improve 3d mesh estimation from video by including audio.
  - E.g. video + audio of ceramic mug being tapped, learn mug mesh, showing improvement over video-only baseline.
  - Downside: Need to re-solve the eigenvalue problem every time the graph changes.
- Learn a 3d mesh that generates audio based on excitations at vertices, only from audio.
  - Use mesh2audio as a prior.
  - Downside: Again - Need to re-solve the eigenvalue problem every time the graph changes.
    - Maybe we can address this problem specifically w/ GNNs?
- Given a 3d mesh and audio sample, estimate the physical material properties (essentially the Faust params)
  - Using mesh2audio model to generate audio.
  - Don't see how GNNs are obviously useful here, but maybe.
- Learn (mesh, excitation index) from: (video, audio, material)
  - Based on the material label, we can set the material props in the mesh2audio modal synthesis model.
  - Use self-supervised training over a completely synthetic dataset of (mesh, material props, excitation index, generated audio clip) samples, as generated by mesh2faust.

Random thoughs:
- Differential DSP suffers from a lack of convexity, easily getting stuck in local minima.
- Learning a 3d mesh that generates audio based from excitations at vertices may be a prior that is more easily optimized, since we can differentiate over 3D space.


### Inspiration

A collection of inspiring papers, not necessarily related to GNNs, but the vibe I'm looking for:
- Differentiable Rendering and Identification of Impact Sounds
  - Paper: https://openreview.net/forum?id=wVIqlSqKu2D
  - Video: https://youtu.be/20-b76wtZ-Q?si=6HElAaDHiC1TZ0MI
- Learning Mesh-Based Simulation with Graph Networks
  - Paper: https://arxiv.org/abs/2010.03409
  - Video: 