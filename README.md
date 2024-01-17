![NMF in Stata](assets/package.png?raw=true "NMF in Stata")

![StataMin](https://img.shields.io/badge/stata-2024-blue) ![issues](https://img.shields.io/github/issues/jonathanbatty/stata-nmf) ![license](https://img.shields.io/github/license/jonathanbatty/stata-nmf) ![version](https://img.shields.io/github/v/release/jonathanbatty/stata-nmf) ![release](https://img.shields.io/github/release-date/jonathanbatty/stata-nmf) ![Stars](https://img.shields.io/github/stars/jonathanbatty/stata-nmf) 

---

[Installation](#Installation) | [Syntax](#Syntax) | [Examples](#Examples) | [Feedback](#Feedback) | [Change log](#Change-log) | [Roadmap](#Roadmap)

---

# Non-negative Matrix Factorisation (NMF) for Stata
(v1.00, 16 Jan 2024)

This repository contains a number of implementations of non-negative matrix factorisation (NMF) using simple multiplicative update rules, implemented in Stata using optimised Mata functions.

## Installation
The package can be installed from GitHub using `net install`:

```
net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

```

## Syntax
The syntax for `nmf` is as follows:

```
nmf varlist, k() [options]

[options] = epoch() initial() loss() stop() nograph()
```

See the help file using `help nmf` for full details of each option.

The most basic usage is as follows:

```
nmf colvar_*, k(n)
```

Whereby a factorisation of the matrix stored in colvar_1, colvar_2, colvar_3, ... , colvar_n will be performed, resulting in matrices W and H of rank n. The resulting matrices will be stored in frames W and H and can be accessed using `frame change W` and `frame change H`, respectively. A summary of the loss parameters by epoch can be viewed using `frame change error`. Note: return to the original dataframe using `frame change default`.

## Examples
Examples of running NMF in 3 different contexts are given in `./examples/`. These use the sample datasets included with the nmf package: `faces.dta`, `nsclc.dta` and `trajectories.dta`.

## Feedback
Please [open an issue](https://github.com/jonathanbatty/stata-nmf/issues) to report errors, suggest feature enhancements, and/or make any other requests. 

## Change Log
**v1.00 (16/01/24)**
 - Initial release.

## Roadmap
- Tentative plan to implement core NMF solver using C++ plugin for speed and stability.

## Acknowledgements
This project was funded by the AI for Science and Government Fund, via the Alan Turing Institute Health Programme (TU/ASG/R-SPEH-114). JB received funding from the Wellcome Trust 4ward North Clinical Research Training Fellowship (227498/Z/23/Z; R127002). 

## Sugested Citation
Batty, J. A. (2024). Stata package ``nmf'': an implementation of non-negative matrix factorisation in Stata (Version 1.0) [Computer software]. https://github.com/jonathanbatty/stata-nmf