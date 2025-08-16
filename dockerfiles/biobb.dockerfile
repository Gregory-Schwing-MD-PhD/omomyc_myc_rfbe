# Base: NVIDIA GROMACS GPU container
FROM nvcr.io/hpc/gromacs:2023.2

# Switch to root to install system tools
USER root

# Install minimal utilities required for Mambaforge
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Download and install Miniforge (Linux installer)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-`arch`.sh -O /tmp/Miniforge3.sh \
    && bash /tmp/Miniforge3.sh -b -p /opt/miniforge \
    && rm /tmp/Miniforge3.sh
ENV PATH=/opt/miniforge/bin:$PATH

# Copy your biobb environment YAML
COPY --chown=root:root envs/biobb_pmx.yml env.yml

# Install environment from YAML
RUN mamba env update --file env.yml -n base \
    && mamba clean --all \
    && mamba update -c conda-forge procps-ng -y \
    && mamba clean --all --yes

# Make the conda environment active by default
SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]

# Default command
CMD ["bash"]
