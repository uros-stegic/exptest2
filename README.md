# exptest2

exptest2 is R package that provides implementations of new (integral- and Kolmogorov- type) exponentiality tests
introduced in paper <sup>[needs citation here]</sup> and also reimplementations of those provided by the package [exptest](https://github.com/cran/exptest).

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

In ordrer to run and develop this software, you'll need the following things installed on your computer:

 - R
 - RStudio (this tutorial uses RStudio, console-type usage is not covered here)
 - rcpp (available in CRAN repository)
 - boost
 - C++ compiler (with C++11 support)
 - rbenchmark (optional, used to measure running times, also available in CRAN)

### Installing

#### Building from source
To build this package, open it in RStudio and press the shortcut Ctrl+Shift+B which will compile the project and open new environment with this library included and ready to use.

#### Using precompiled binaries
First, you need to download the binary package [here](https://github.com/uros-stegic/exptest2/releases/tag/v0.1.2-alpha). After that, open RStudio and from the menu bar select `Tools -> Install Packages...` and in the popup window, in the `Install From:` dropdown menu select `Package Archive File (.tar.gz)`. Then find and select the binary package you downloaded and hit `Install`.

## Usage
ExpTest exports two functions for testing wether the given distribution comes from exponential family. See this example:
```R
# Shows usage of both exp_test_integral and exp_test_kolmogorov functions
> library(ExpTest)
> x <- rexp(50)
> exp_test_integral(x)
[1] 0.003244134
> exp_test_kolmogorov(x)
[1] 0.07005984
```
For speed testing, you can use `rbenchmark` library that can be installed from CRAN. To measure speed of these functions, simply call the `benchmark` function like this:
```R
# Measure running times of these two functions
> library(ExpTest)
> library(rbenchmark)
> x <- rexp(30)
> benchmark(replications = rep(1, 1, 1), exp_test_integral(x), exp_test_kolmogorov(x))
                    test replications elapsed relative user.self sys.self user.child sys.child
1   exp_test_integral(x)            1   0.018        1     0.019        0          0         0
2 exp_test_kolmogorov(x)            1   0.018        1     0.018        0          0         0
```


## Authors

* **Uros Stegic** - *Initial work* - [uros-stegic](https://github.com/uros-stegic)

## License

This project is licensed under the GNU GPL-3 License - see the [LICENSE](LICENSE) file for details

