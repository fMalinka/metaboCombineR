# metaboCombineR - A tool for combining untargeted multi-batch LC-MS experiments

metaboCombineR is R package written in C++ that focuses on large-scale untargeted LC-MS experiments. The package allowes to preprocessed each batch of samples separately using an arbitrary preprocessing program, such as XCMS. The main pros is a possibility to handle and observed results of experiments during the time and program parameter tuning should be easier since within-batch variablity is smaller than beetween-batch variablity.

## Getting Started

### Prerequisites
Only one external R package is required: `Rcpp (>= 0.12.16)`. For compiling vignette we suggest `rmarkdown` package.
### Installing
The simplest way to install metaboCombineR package is via `devtools` that allowes to download and install the project instantly from gitHub using only one command.
```
> library(devtools)
> install_github("fmalinka/metaboCombineR")
```

## Running the example
For an ilustration, we prepared four authentic real datasets. For loading them to workspace, type the following:
```
> library(metaboCombineR)
> data(metaboExp1)
> data(metaboExp2)
> data(metaboExp3)
> data(metaboExp4)
```
Presented datasets are in 2-dimensional matrix format where rows represent features and each row has it own name which m/z value is prefixed by `M` and rt by `T`. Columns represent samples. Names of all samples must be filled by `colnames` function through R. Values of matrix must be in numerical format, not `factor` type.

To combine all of these experiments into one table, call `runMetaboCombiner` function, where the first arguments is supposed to be a list of experiments, `mzprecision` argument defines a number of digits considered for peaks.
```
mytable <- runMetaboCombiner(list(metaboExp1, metaboExp2, metaboExp3, metaboExp4), mzprecision = 2)
```
Then, the result matrix is stored in `mytable` variable.

## Authors

* **František Malinka**

## Citation
Cite please:


Malinka, F., Zareie, A., Novosadova, V. Batch alignment via retention orders for improving the reliability of untargeted metabolomics. Unpublished.

## Acknowledgments

* Ashkan Zareie for his considered comments and suggestions.

## License

This project is licensed under the MIT License.

