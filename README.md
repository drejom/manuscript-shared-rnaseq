RNA-Seq code and resources used in multiple manuscripts
=======================================================

## Overview

This repository contains a collection of R code, [RMarkdown documents](http://rmarkdown.rstudio.com/),
and data/annotations which have been used for several different Host-pathogen 
RNA-Seq analyses performed by our lab.

The main components of this repository are organized as follows:

- `data/` - Data and annotations generated in our lab, or accessed from
  open-access journal articles of interest.
- `R/` - R functions for differential expression analysis, plotting, etc. which
  are generally useful across multiple analyses.
- `Rmd/` - RMarkdown documents each performing some specific task, e.g.
  RNA-Seq count filtering.

The `data/` section is further sub-divided by species for which the resources
are relevant, while the `R/` and `Rmd/` folders are divided by the function or
task to which they relate.

The RMarkdown documents have been designed in such a way so that each document
performs a single function, or some small set of functions, and can be imported
as [knitr child documents](http://yihui.name/knitr/demo/child/) into a larger
manuscript-specific analysis.

Individual projects/manuscripts can then build on top of these modular `.Rmd`
child documents, picking-and-choosing the components which are relevant to
those analyses and extending them with experiment-specific functionality.

## Functionality

The code and functionality within this repository is largely geared towards two
types of analysis:

- Differential Expression Analysis
- Co-expression Network Analysis

In addition to these two broad types of analysis, the repo also includes code
for handling:

- Batch adjustment and normalization
- Visualization of RNA-Seq sample relationships
- RNA-Seq count filtering
- Enrichment analysis (GO, KEGG, ConsensusPathDB, LeishCyc, ENCODE)

## Usage

### Overview

To make use of this functionality you must:

1. Create a parent Rmarkdown document
2. Clone this repository to an accessible location
2. Define a `CONFIG` list variable containing the analysis settings
3. Use `knit_child()` or the `child=` chunk option to load individual Rmarkdown
   child documents in the desired order.
5. Use Rstudio or the `render()` function in the [rmarkdown
   package](https://cran.r-project.org/web/packages/rmarkdown/index.html) to
   render this document into HTML, PDF, etc.

### Example analysis

For example, to perform a basic differential expression analysis one might
create a parent `manuscript.Rmd` document which contains:

    ---
    title: "Example Differential Expression Manuscript"
    ---

    ```{r child='shared/Rmd/init/header_de.Rmd'}
    ```

    ```{r child='analysis_settings_v1.0.Rmd'}
    ```


    ```{r child='shared/Rmd/init/load_settings_de.Rmd'}
    ```

    Methods
    =======

    ```{r child='shared/Rmd/init/load_counts.Rmd'}
    ```

    ```{r child='shared/Rmd/init/load_host_annotations.Rmd'}
    ```

    ```{r child='shared/Rmd/init/create_expression_set.Rmd'}
    ```

    ```{r child='shared/Rmd/main/filter_counts.Rmd'}
    ```

    ```{r child='shared/Rmd/main/gene_visualization.Rmd'}
    ```

    ```{r child='shared/Rmd/main/data_prep_de.Rmd'}
    ```

    ```{r child='shared/Rmd/main/differential_expression.Rmd'}
    ```

    ```{r child='shared/Rmd/results/go_enrichment_de.Rmd'}
    ```

    ```{r child='shared/Rmd/results/kegg_enrichment_de.Rmd'}
    ```

    ```{r child='shared/Rmd/results/save_output_de.Rmd'}
    ```

    Version Information
    ===================

    ```{r}
    sessionInfo()
    ```

### Configuration

In the above example, all but one of the child documents which get imported are
contained within this repository.

The one file which is not part of this repo is the file called "analysis_settings_v1.0.Rmd".

This is where all of your experiment-specific settings and desired analysis
parameters can be specified.

This includes things like:

- the name of the host and pathogen species being analyzed
- data pre-processing and normalization settings
- p-value cutoffs
- etc.

All of these settings are specified in a single global list variable named `CONFIG`,
for example:

    ```{r}
    config <- list(
        # general
        'pathogen'            = 'l. major',
        'host'                = 'h. sapiens',
        'target'              = 'host',
        'analysis_name'       = 'hsapiens-inf-with-lmajor-1.0',

        # sample variables of interest
        'sample_id'           = 'sample_id',
        'condition'           = 'condition',
        'batch'               = 'replicate',

        # sample metadata
        'samples'             = tbl_df(read.csv('hsapiens_samples.csv'))

        # input count tables
        "input_count_tables"  = file.path("counts/mmusculus/*.count"),

        # normalization / data transformation
        'de_cpm'                 = true,
        'de_log2'                = true,
        'de_voom'                = true,
        'de_quantile_normalize'  = true,
        'de_batch_adjust'        = 'combat',

        # annotations
        "organism_db"         = "homo.sapiens",
        "orgdb_key"           = "ensembl",
        "organism_genome"     = "hg38",
        "biomart_dataset"     = "hsapiens_gene_ensembl",

        # differential expression
        "contrast_column"  = 'condition',
        "de_comparisons"   = list(c('uninfected', 'infected')),
    )
    ```
Where possible, reasonable default options have been chosen so that it is not
necessary to set every single option, but this is not always possible.

For a complete list of settings, you can find the default configuration options
listed in `Rmd/settings/shared.Rmd`, `Rmd/settings/coex_network.Rmd`, and
`Rmd/settings/differential_expression.Rmd`.

You can also refer to the manuscript specific repositories which show actual
example of how this code is being used.

## Bioconductor Docker 

## Building Docker image

To build the image, navigate to root directory with the `Dockerfile` and run

``` sh
time docker build -t denom/bioconductor_docker_genomics:latest . 2>&1 | tee docker.build.log
```

Note errors:

``` sh
grep "ERROR:" docker.build.log
```

## Pushing the image to DockerHub

First login to Docker on the command line. Find the hash of the image

`cat ~/.docker/my_password.txt | docker login --username denom --password-stdin`

Then push to DockerHub

```sh
docker push denom/bioconductor_docker_genomics:latest
```

## Running the container

```sh
docker run -d \
   -v $(pwd):/home/rstudio/data \
   -v $HOME/R:/usr/local/lib/R/host-site-library \
   -e PASSWORD=bioc \
   -p 6080:8787 \
   bioconductor/bioconductor_docker_genomics:RELEASE_3_10
```

Contact Information
===================

For questions or suggestions related to this analysis, you are welcome to other
submit fixes to Github directly, or contact me at
[khughitt@umd.edu](mailto:khughitt@umd.edu).


