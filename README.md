# Multigrid OOP Solver
1D/2D/3D finite difference multigrid solver on regular grid. The skeleton of the code is the same as the perfect 2D multigrid solver provided by Achi Brandt. Everyone who is new to this should learn how masters do their work. I personally definitely want to meet the author and learn from him.

This code also serves as an example of how one can simply/generalize an existing code. Original documents can be downloaded [here](http://www.wisdom.weizmann.ac.il/~achi/classics.pdf). Useful references are: [Matlab code tutorial](https://github.com/dappelha/MultiGridMatlab)

## Getting Started

This code is completely in OOP. Choose the dimension you are aiming for, and start from:
```
TestCycle.run;
```

All the parameters are set in 'Options.m'. Currently weight Jacobi and two types of Gauss-Seidal smoother are implemented.

### Prerequisites

MATLAB

## Notes
There are rooms for improvements. I found a [master thesis](https://www.duo.uio.no/handle/10852/12685) on this topic which a good summary and review of the related methods. One key idea from that thesis is to use convolution instead of explicit loops over indexes to do relaxation. Based on my tests, this gives faster speed for a grid size larger than 64, but it is not worth it for any smaller grid. The current 3D implementation of GS method is a combination of the brute-force method and convolution method.

Line restrictor and interpolator are implemented. This works especially well if there is a strong decoupling between different directions. For 3D, it is also possible to use plane smoother.

Alternatively, other 3D matrix stencil restrictor and interpolator are implemented and tested. I have no idea which is the best --- the usual answer to these kind of questions is: it depends.

Currently only Dirchlet boundary condition is applied. It shouldn't be difficult to add other common boundary conditions.

## Authors

* **Hongyang Zhou** - [henry2004y](https://github.com/henry2004y)
* **Achi Brandt** - *Original code*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* All the nice guys who share their codes


