# blockRAR: Block Design for Response Adaptive Randomization

[![Build Status](https://travis-ci.org/thevaachandereng/blockRAR.svg?branch=master)](https://travis-ci.org/thevaachandereng/blockRAR)
[![Download badge](https://cranlogs.r-pkg.org/badges/blockRAR)](https://cran.r-project.org/package=blockRAR)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![codecov](https://codecov.io/gh/thevaachandereng/blockRAR/branch/master/graph/badge.svg)](https://codecov.io/gh/thevaachandereng/blockRAR)

**Authors**: Thevaa Chandereng and Rick Chappell


Overview
--------
Response-Adaptive Randomization (RAR) is an adaptive trial where the randomization ratio of the patient changes based on the performance of the control and experimental treatment. 
However, most design complete ignores the time trend aspect in this design and the randomization
ratio are altered based on a patient outcome. 
blockRAR assigns patient in a block (group) manner and the the block results are analyzed before the randomization ratio is altered.
Time is divided into factor level in each block (group).
The treatment and time effect is both obtained in this design. 
The blockRAR website is available [here](https://thevaachandereng/blockRAR/). 


Installation
------------
Prior to analyzing your data, the R package needs to be installed.

The easiest way to install blockRAR is through CRAN:

``` r
install.packages("blockRAR")
```

There are other additional ways to download blockRAR.
The first option is most useful if want to download a specific version of blockRAR
(which can be found at https://github.com/thevaachandereng/blockRAR/releases).
``` r 
devtools::install_github("thevaachandereng/blockRAR@vx.xx.x")
# OR 
devtools::install_version("LPWC", version = "x.x.x", repos = "http://cran.us.r-project.org")
```

The second option is to download through GitHub. 

``` r
devtools::install_github("thevaachandereng/blockRAR")
```

After successful installation, the package must be loaded into the working space:

``` r 
library(blockRAR)
```

Usage
------------
See the [vignette](https://gitter-lab.github.io/LPWC/articles/LPWC.html) for usage instructions.


Reference
------------
If you use blockRAR, please cite TBA. 

License
------------
blockRAR is available under the open source [MIT license](http://opensource.org/licenses/MIT).