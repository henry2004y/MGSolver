# Multigrid OOP Solver
1D/2D/3D finite difference multigrid solver on regular grid. The skeleton of the code is the same as the perfect 2D multigrid solver provided by Achi Brandt. Everyone who is new to this should learn how masters do their work. I personally definitely want to meet the author and learn from him.

This code also serves as an example of how one can simply/generalize an existing code. Original documents can be downloaded [here](http://www.wisdom.weizmann.ac.il/~achi/classics.pdf).

## Getting Started

This code is completely in OOP. Choose the dimension you are aimed for, and start from:
```
TestCycle.run;
```

All the parameters are set in 'Options.m'. Currently weight Jacobi and two types of Gauss-Seidal smoother are implemented.

### Prerequisites

MATLAB

## Notes
There are rooms for improvements. I found a master thesis on this topic which is very useful.

## Authors

* **Hongyang Zhou** - [henry2004y](https://github.com/henry2004y)
* **Achi Brandt** - *Original code*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* All the nice guys who share their codes


