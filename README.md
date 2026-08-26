t# phybreak
Outbreak reconstruction with sequence data

The package implements the method described in Van der Roest et al. (2023), https://doi.org/10.1371/journal.pcbi.1010928
Workflow:

* enter data and priors by constructing an object of S3-class 'phybreak', with function 'phybreak'

* do mcmc-updates with functions 'burnin_phybreak' and 'sample_phybreak'; remove samples with 'thin.phybreak'

* access the 'phybreak'-object by get_phybreak-functions such as 'get_transtree', 'get_data', 'get_parameters'

* summarize the mcmc-chain with the functions 'ESS', 'transtree', 'infectorsets', 'phylotree', 'get_mcmc', 'get_phylo'

* plotting with 'plot', 'plotTrans', and 'plotPhylo'


* it is possible to simulate data with 'sim_phybreak'

## Installation

The package requires some external libraries.
Please install them by typing:

```bash
sudo apt-get install libblas-dev liblapack-dev
```

The package can be directly installed from GitHub using the `devtools` R package.

```r
devtools::install_github("https://github.com/bastiaanvdroest/phybreak/")
```

## Contact type contributions

To estimate the contributions of different contact types on transmission, use the code on the contact branch of this repository.
Here, the function 'phybreak' has the option to select 'contact = TRUE', such that the contact contributions are taken into account in the likelihood calculations of the MCMC algorithm in the model.
