# Modeling & Simulation Project

## Plan

* Ben:
  - Render an input volumetric mesh to mass/stiffness eigenvectors & eigenvalues.
  - Everything in the "Vega FEM" box (in the `mesh2faust` implementation diagram below) needed to implement this specific FEM model:
    * Given volumetric mesh (in [GLTF 2.0 format](https://github.com/KhronosGroup/glTF)) & material properties, produce:
      - Tetrahedral mesh -> (Mass matrix, Stiffness matrix) -> Eigen solver -> (Eigenvectors, Eigenvalues)
* Karl:
  - Render the FEM outputs to a modal physical audio model.
  - Everything outside of the "Vega FEM" box:
    * Provide volumetric mesh & material properties to FEM.
    * UI/UX/audio/mesh generation & editing
    * Take eigenvectors and eigenvalues from "Vega FEM" output, translating them to modes gain/frequency matrices.
    * Audio model parameter editor
    * Given model parameters + mass/stiffness eigenvectors/eigenvalues:
      - Produce (Gain/frequency matrices) -> Modes selection -> Audio rendering
    * Continuously fill the audio buffer based on the current mesh modes.
      - Initially implement using `mesh2faust->faust`, JIT rendering Faust to LLVM for runtime audio graph updates.
      - Next, render Faust->Julia->LLVM, JIT-loading LLVM on changes.
      - Next, render mesh modes directly to Julia, without Faust inbetween.

## Stack

* FEM: Finite element model of 
* App: UI/UX/audio/mesh generation & editing
  * [ImGui](https://github.com/ocornut/imgui): Renders the UI to many backends.
  * [GLTF 2.0](https://github.com/KhronosGroup/glTF): Khronos spec. Transmits and loads 3D scenes and models efficiently.
    Minimizes the size of 3D assets, and their unpack runtime.
  * [Filament](https://github.com/google/filament): Real-time physically based rendering library.
    GLTF loader/viewer/editor.
  * [MeshOptimizer](https://github.com/zeux/meshoptimizer): Optimizes meshes for GPU pipelining, and reduces mesh complexity & storage overhead.
  * [Faust](https://github.com/grame-cncm/faust): Renders the mesh to an audio graph, with real-time interactive vertex exitation.
  * [miniaudio](https://github.com/mackron/miniaudio): Continuously renders the modal physical model of the input 3D volumetric mesh to audio.

## M&S project direction

It would be better to keep the project standalone.
- The whole damn thing is in one place, and my project isn't seen or beholden to anything.
- I get more practice building apps in this stack from scratch.
- More experience treating these libraries like most users do, without the wrapper semantics of FlowGrid.
- More immediacy. Can't possibly expect a teammate to learn about FlowGrid. Ok that's it, done.
- More fun!!!
- Can (and will!) use it in FlowGrid anyway!
  * Suggests designing like a library.

## Development strategy

* Build it like a library.
* Make it performant.
* Use a single `struct State` intance
* Communicate by updating state.
* Prefer read-only members to getter methods
* Don't even think about working on app state management.
  - State import/export only
  - (State management remains a FlowGrid focus.
    It can use this as a mesh editor library later.)
