FROM rocker/r-ver:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libglpk40 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    libreadline-dev \
    libxt-dev \
    pandoc \
    && apt-get clean

# Install CRAN packages
RUN Rscript -e "install.packages(c( \
  'Seurat', \
  'ggplot2', \
  'SeuratObject', \
  'remotes', \
  'devtools', \
  'rmarkdown', \
  'knitr', \
  'Matrix', \
  'dplyr', \
  'data.table' \
), repos = 'https://cloud.r-project.org')"

# Install Bioconductor packages
RUN Rscript -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')" \
    && Rscript -e "BiocManager::install(c( \
      'SingleR', \
      'celldex', \
      'SummarizedExperiment', \
      'rtracklayer' \
    ), ask = FALSE, update = FALSE)"

# Copy the local package into the container
COPY . /home/rstudio/programmingpackage

# Install the custom package
RUN Rscript -e "remotes::install_local('/home/rstudio/programmingpackage', force = TRUE, upgrade = 'never')"

# Set the working directory
WORKDIR /home/rstudio

# Start R by default
CMD ["R"]

