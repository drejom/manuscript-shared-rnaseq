FROM bioconductor/bioconductor_docker:RELEASE_3_10

LABEL name="denom/bioconductor_docker_genomics" \
      version="0.0.1" \
      url="https://github.com/drejom/bioconductor_docker_genomics" \
      maintainer="omeally@gmail.com" \
      description="Description of my image" \
      license="Artistic-2.0"

# Install dependencies
RUN apt-get update && apt-get install -y cargo \
    && apt-get clean \                 
    && rm -rf /var/lib/apt/lists/*


## Install R packages 
ADD install.R /tmp/
RUN R -f /tmp/install.R

# Open a port for RStudio server
EXPOSE 8787
