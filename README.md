# Shortest-Path Percolation on Random Networks

This repository contains the code for
- Minsuk Kim and Filippo Radicchi, [Shortest-Path Percolation on Random Networks](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.133.047402), Physical Review Letters $133$, 047402 (2024).
- [Preprint (arXiv)](https://arxiv.org/abs/2402.06753)

- BibTex entry:
    ```
    @article{kim2024shortest,
    title={Shortest-Path Percolation on Random Networks},
    author={Kim, Minsuk and Radicchi, Filippo},
    journal={Physical Review Letters},
    volume={133},
    number={4},
    pages={047402},
    year={2024},
    publisher={APS}
    }
    ```


# How to Run the Code
We provide a `makefile` that compiles the source code. You can compile it by typing:

```
make
```

Then, you can run the code as follows:
```
efficient_SPP_ERN.out <N> <avg_k> <C> <num_iter> <num_instance>
```

where `N` is the size of the ER graph, `avg_k` is the mean degree of the ER graph, `C` is the maximum shortest-path length, `num_iter` is the number of iterations of the shortest-path percolation on an instance of the ER graph, and `num_instance` is the number of ER graph instances.

We also provide a Jupyter notebook file, `example.ipynb`, to demonstrate some basic analysis of the output data from `efficient_SPP_ERN.out`.
