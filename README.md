![NMF in Stata](assets/package.png?raw=true "NMF in Stata")

![StataMin](https://img.shields.io/badge/stata-2024-blue) ![issues](https://img.shields.io/github/issues/jonathanbatty/stata-nmf) ![license](https://img.shields.io/github/license/jonathanbatty/stata-nmf) ![Stars](https://img.shields.io/github/stars/jonathanbatty/stata-nmf) ![version](https://img.shields.io/github/v/release/jonathanbatty/stata-nmf) ![release](https://img.shields.io/github/release-date/jonathanbatty/stata-nmf)

---

[Installation](#Installation) | [Syntax](#Syntax) | [Examples](#Examples) | [Feedback](#Feedback) | [Change log](#Change-log) | [Roadmap](#Roadmap)

---

# Non-negative Matrix Factorisation (NMF) for Stata, v1.00
(16 Jan 2024)

This repository contains a number of implementations of non-negative matrix factorisation (NMF) using simple multiplicative update rules, implemented in Stata using optimised Mata functions.

## Installation
The package can be installed from GitHub, using the `net install` syntax:

```
net install nmf, from("https://raw.githubusercontent.com/jonathanbatty/stata-nmf/main/installation/") replace

```

## Syntax
Syntax

## Examples
Examples

## Feedback
Please open an [issue](https://github.com/jonathanbatty/stata-nmf/issues) to report errors, feature enhancements, and/or other requests. 

## Change Log
**v1.00 (16/01/24)**
 - Initial release.


## Roadmap
- Re-implement solver using C++ for speed and stability.

## Acknowledgements
This project was funded by the AI for Science and Government Fund, via the Alan Turing Institute Health Programme (TU/ASG/R-SPEH-114). JB received funding from the Wellcome Trust 4ward North Clinical Research Training Fellowship (227498/Z/23/Z; R127002). 

## Sugested Citation
Batty, J. A. (2024). Stata package ``nmf'': an implementation of non-negative matrix factorisation in Stata (Version 1.0) [Computer software]. https://github.com/jonathanbatty/stata-nmf