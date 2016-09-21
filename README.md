# ControlCore

Core functionality for **identifying**, **analyzing** and **designing** control
systems in `Julia`.

### Build Status and Code Coverate

-  Build status: [![Build Status][build-img]][build-link]
-  Code coverage: [![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

[build-img]:  https://travis-ci.org/KTH-AC/ControlCore.jl.svg?branch=master
[build-link]: https://travis-ci.org/KTH-AC/ControlCore.jl
[ca-img]: https://coveralls.io/repos/github/KTH-AC/ControlCore.jl/badge.svg?branch=master
[ca-link]: https://coveralls.io/github/KTH-AC/ControlCore.jl?branch=master
[cc-img]: https://codecov.io/gh/KTH-AC/ControlCore.jl/branch/master/graph/badge.svg
[cc-link]: https://codecov.io/gh/KTH-AC/ControlCore.jl

### Description

This repository is meant to provide a basic set of tools, *i.e.*, *transfer
function* and *state space* types as well as the basic mathematical operations
defined on them, in a way that other tools for [identifying][kth-sysid] or [analyzing][kth-ct] control systems would use the same set of functionality
in a transparent, coherent way.

[kth-sysid]: https://github.com/KTH-AC/IdentificationToolbox.jl
[kth-ct]: https://github.com/KTH-AC/ControlToolbox.jl
