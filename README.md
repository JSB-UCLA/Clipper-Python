[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg )](https://github.com/JSB-UCLA/Clipper-Python/blob/main/LICENSE.md)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

# Clipper-Python

The Python package for Clipper: A p-value-free method for controlling false discovery rates in high-throughput biological data with two conditions.

**Clipper** aims to contrast two conditions to reliably screen interesting features. The interesting means "differential" or "enriched". Differential
features are defined as those that have different expected measurements between two conditions, and enriched features are those that have higher expected measurements under the experimental/treatment condition than the background condition, i.e., the negative control. 

Any suggestions on the package are welcome! For suggestions and comments on the method, please contact Xinzhou (<xinzhouge@ucla.edu>) or Dr. Jessica Li (<jli@stat.ucla.edu>).

* **Free software:** MIT license
* **Python versions:** 3.6 and above

## Installation

To install Clipper-Python, clone this repo:

```shell
$ git clone https://github.com/JSB-UCLA/Clipper-Python.git
```
and then run

```shell
pip install -r ./Clipper/requirements.txt
```
to install required packages

To import Clipper, you need to first add the directory of the Clipper module into the system path:

```python
import sys
base_dir = "dir/to/Clipper"
sys.path.append(base_dir)
import Clipper
```


### All-in-one function

Clipper requires a minimum of four inputs:

- `score_exp` a numeric matrix of measurements under the experimental condition with rows being features and columns being replicates.
- `score_back` a numeric matrix of measurements under the background condition with rows being features and columns being replicates. Its features should match features from `score_back`.
- `analysis` a character string specifying the analysis goal, must be either "differential" ("d") or "enrichment" ("e"). See below for details.
- `FDR` a numeric value or vector indicating the target FDR threshold(s).

Optional input arguments include:

- `procedure` a character string specifying the FDR control procedure, must be "GZ", "BC" or "aBH". If not specified, Clipper uses the GZ procedure when `analysis = "differential"` and the BC procedure when `analysis = "enrichment"`
- `contrast.score` a character string specifying the contrast score, must be "max" or "diff". If not specified,
`Clipper` uses the maximum (max) contrast score when `procedure = "GZ"` and the difference (diff) contrast score when `procedure = "BC"`.
- `n.permutation` an integer specifying the number of permutations. Effective only when `procedure = "GZ"`
- `seed` random seed, used in permutations.

Clipper returns a list containing the following components:
- `contrast.score` the type of contrast scores, 'max' or 'diff';
- `contrast.score.value` the values of contrast scores;
- `FDR` a vector of the target FDR threshold(s);
- `contrast.score.thre` a vector of the threshold(s) on contrast scores corresponding to the FDR threshold(s);
- `q` the q-values of each feature;
- `discovery` a list of identified discoveries. Each component contains discovered features, coded as the row indices of `score_exp` and `score_back`, at a FDR threshold.

You can use the following command to run clipper function once the Clipper module is imported:
```python
re = Clipper.clipper(score_exp=exp_d, score_back=back_d, analysis="differential", FDR=[0.05])
```

### Example

You can find examples of Clipper usage in different settings in "test_Clipper.py".
