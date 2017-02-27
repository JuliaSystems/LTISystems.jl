# SystemsBase

Core functionality for **identifying**, **analyzing** and **designing** control
systems in `Julia`.

### Build Status and Code Coverage

-  Build status: [![Build Status][build-img]][build-link]
-  Code coverage: [![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

[build-img]:  https://travis-ci.org/JuliaSystems/SystemsBase.jl.svg?branch=master
[build-link]: https://travis-ci.org/JuliaSystems/SystemsBase.jl
[ca-img]: https://coveralls.io/repos/github/JuliaSystems/SystemsBase.jl/badge.svg?branch=master
[ca-link]: https://coveralls.io/github/JuliaSystems/SystemsBase.jl?branch=master
[cc-img]: https://codecov.io/gh/JuliaSystems/SystemsBase.jl/branch/master/graph/badge.svg
[cc-link]: https://codecov.io/gh/JuliaSystems/SystemsBase.jl

### Description

This repository is meant to provide a basic set of tools, *i.e.*, *transfer
function* and *state space* types as well as the basic mathematical operations
defined on them, in a way that other tools for [identifying][kth-sysid] or [analyzing][kth-ct] control systems would use the same set of functionality
in a transparent, coherent way.

[kth-sysid]: https://github.com/JuliaSystems/IdentificationToolbox.jl
[kth-ct]: https://github.com/JuliaSystems/ControlToolbox.jl
