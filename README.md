[![DOI](https://zenodo.org/badge/151991844.svg)](https://doi.org/10.5281/zenodo.15033018)

# R Functions for Hydrology and Morphometry

This repository contains a collection of R functions specifically designed for computational hydrology and morphometric analyses. These functions leverage the capabilities of GRASS GIS 7 and utilize the archived R package `rgrass7`. They are particularly suited for detailed studies of watersheds and drainage networks.

## Features

- **Hydrological Analysis**: Tools for analyzing water flow, watershed delineation, and drainage patterns.
- **Morphometric Analysis**: Functions to assess terrain shapes, elevation profiles, and other geomorphological features.
- **Integration with GRASS GIS 7**: Seamless use of GRASS GIS functionalities through R.

## Getting Started

### Prerequisites

Before you can use these functions, you need to have GRASS GIS 7 installed on your computer. Additionally, the R package `rgrass7` must be installed and configured to interact with your GRASS GIS installation.

### Installation

Download this repository to your local machine using the `Code>Download ZIP` button of the repo.

Navigate to the downloaded ZIP and extract it.

Install the required R package `rgrass7` by downloading it from the [archived version on CRAN](https://cran.r-project.org/web/packages/rgrass7/index.html)

## Usage

To use the provided functions, load them into your R session:

```R
source("path/to/function.R")
```

For example:

```R
source("rgrass/lfp_network.R")
```

Examples and detailed usage instructions are provided in the documentation within each function file.

## Acknowledgments

- Thanks to the GRASS GIS team for their powerful GIS platform.
- This project builds upon the work done in the `rgrass7` package, providing additional functionality tailored for hydrological and morphometric analysis.
