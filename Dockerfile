# Use the litd/docker-cellranger as the base image
FROM litd/docker-cellranger

# Set the working directory and run the following commands there
WORKDIR /home/loka
# Set an environment variable for non-interactive installs
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN dnf -y update && \
    dnf -y install \
    gcc \
    gcc-c++ \
    make \
    wget \
    git \
    curl \
    openssl-devel \
    libffi-devel \
    python3-devel \
    && dnf clean all

# Install Miniconda (used to install conda and mamba)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set the path for conda
ENV PATH="/opt/conda/bin:$PATH"

# Install Mamba into the base environment
RUN conda install -n base -c conda-forge mamba

# Copy the condaenv.yml file to the container
COPY condaenv.yml /tmp/condaenv.yml

# Use Mamba to create the environment from the .yml file
RUN mamba env create -f /tmp/condaenv.yml

# Clean up after installation to reduce image size
RUN dnf clean all && \
    conda clean --all --yes

# Activate the environment by default
# Make sure the environment name matches the one in your .yml file
ENV PATH="/opt/conda/envs/condaenv/bin:$PATH"
SHELL ["conda", "run", "-n", "condaenv", "/bin/bash", "-c"]

# Set the default command to bash
CMD ["/bin/bash"]
