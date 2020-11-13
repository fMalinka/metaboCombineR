# metaboCombineR - A tool for combining untargeted multi-batch LC-MS experiments

metaboCombineR is R package written in C++ that focuses on large-scale untargeted LC-MS experiments. The package allowes to preprocessed each batch of samples separately using an arbitrary preprocessing program, such as XCMS. The main pros is a possibility to handle and observed results of experiments during the time and program parameter tuning should be easier since within-batch variablity is smaller than beetween-batch variablity. The packages contains two different approaches called `kmersAlignment` and `rtcorrectedAlignment` that were both published in [1].

## Getting Started

### Prerequisites
Only one external R package is required: `Rcpp (>= 0.12.16)`. For compiling vignette we suggest `rmarkdown` package.
### Installing
The simplest way to install metaboCombineR package is via `devtools` that allowes to download and install the project instantly from gitHub using only one command.
```
library(devtools)
install_github("fmalinka/metaboCombineR")
```

## Running the example
For an ilustration, we prepared four authentic real datasets. For loading them to workspace, type the following:
```
library(metaboCombineR)
data(metaboExp1)
data(metaboExp2)
data(metaboExp3)
data(metaboExp4)
```
Presented datasets are in 2-dimensional matrix format where rows represent features and each row has it own name which m/z value is prefixed by `M` and rt by `T`. Columns represent samples. Names of all samples/features must be filled by `colnames`/`rownames` function through R. To see an example of dataset format, type `head(metaboExp1)`.

MetaboCombineR enables to define assignments of samples to some group, typically as treatment and control groups. To provide this information, add a vector of assignments as the first row to the input `data.frame`. The row must be named as `group`. If the group label is present in a dataset it must be present in all other datasets. Inconsistencies are not allowed. An example of dataset with three samples (two treatments and one control) and 2 features is depicted below.

```
| X064.EPK83_m_Mzb1_ESI.mzML | X064.EPK88_m_Mzb1_ESI.mzML | X064.EPK94_m_Mzb1_ESI.mzML |
| group | treatment | treatment | control |
| M57.08131T1428.18786 | -0.321542424867644 | 0.286250559905367 | 1.17078411764221 |
| M57.23559T1428.09065 | 0.408813652123444 | -1.0100177456997 | -0.153473421445439 |
```

To combine all of these experiments into one table, call `runMetaboCombiner` function, where the first arguments is supposed to be a list of experiments, `mzprecision` argument defines a number of digits considered for peaks, and `algorithm` selects one of the proposed algorithms. The `algorithm` argument supposes only two values on input: kmer and rtcor. The default algorithm is kmersAlignment. `windowsize` argument represents the size of the window for `rtcorrectedAlignment` algorithm. For `kmersAlignment` algorithm, the same argument represents the k-mer value.
```
mytableKmer <- runMetaboCombiner(list(metaboExp1, metaboExp2, metaboExp3, metaboExp4), mzprecision = 2, algorithm="kmer", windowsize = 5)
mytableRtcor <- runMetaboCombiner(list(metaboExp1, metaboExp2, metaboExp3, metaboExp4), mzprecision = 2, algorithm="rtcor", windowsize = 50)
```
Then, the result matrix is stored in `mytableKmer`, `mytableRtcor` variable respectively.

## Authors

* **František Malinka**

## Citation
Cite please:


[1] Malinka, F., Zareie, A., Procházka, J., Sedláček, R., Novosadova, V. Batch alignment via retention orders for preprocessing large-scale multi-batch LC-MS experiments. Unpublished.

## Acknowledgments

* Ashkan Zareie for his considered comments and suggestions.

## License

This project is licensed under the MIT License.

