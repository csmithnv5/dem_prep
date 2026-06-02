# Start from a core stack version
FROM ubuntu:24.04

# Avoid interactive prompts during apt install
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies needed for Miniconda and package builds
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    git \
    && rm -rf /var/lib/apt/lists/*


# Install Miniconda to /opt/conda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

# Put conda on the system PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# Copy dependency configuration files into the container
# (Choose the file type that matches your project below)
COPY environment.yml requirements.txt* ./

RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main \
 && conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
 
# Conda environment.yml file
RUN conda env create -f environment.yml -y && \
	conda clean -afy 

WORKDIR /home

# Ensure the conda environment is activated for the terminal session
RUN echo "conda activate grassenv" >> ~/.bashrc

# Set the default shell to use the conda environment's python bin path
ENV PATH /opt/conda/envs/grassenv/bin:$PATH

# Make RUN commands inside the Dockerfile use the new environment automatically
# Replace 'myenv' with the name defined inside your environment.yml
SHELL ["conda", "run", "-n", "grassenv", "/bin/bash", "-c"]

# Expose Jupyter port
EXPOSE 8888

# Start Jupyter Server inside the Conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "grassenv", "jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''"]