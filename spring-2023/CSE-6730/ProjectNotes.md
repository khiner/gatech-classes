# Modeling & Simulation Project

GH: https://github.com/GATech-CSE-6730-Spring-2023-Project/mesh2audio

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


## Mesh viewer

Got a fork of [mesh-viewer](https://github.com/amanshenoy/mesh-viewer) working on my machine here: https://github.com/khiner/mesh-viewer

* Would be nice to decouple from GLFW, and work with sdl3.
* Port things over one by one?