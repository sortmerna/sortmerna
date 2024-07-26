
FROM ubuntu:jammy

ENV CONDA_DIR=/opt/conda
ENV SMR_DIR=/sortmerna
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND noninteractive
ENV PATH=${SMR_DIR}/dist/bin:${CONDA_DIR}/bin:${PATH}

# jammy packs gcc 11.4.0 by default (build-essential)
RUN apt-get update > /dev/null && \
    apt-get install --no-install-recommends --yes \
        build-essential wget bzip2 ca-certificates git tini file > /dev/null && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget --no-hsts --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh && \
    mamba install -y pyyaml jinja2 requests ninja cmake=3.29.6  && \
    git clone https://github.com/biocore/sortmerna.git && \
    cd ${SMR_DIR} && \
    python setup.py -n all && \
    conda clean --tarballs --index-cache --packages --yes && \
    find ${CONDA_DIR} -follow -type f -name '*.a' -delete && \
    find ${CONDA_DIR} -follow -type f -name '*.pyc' -delete && \
    conda clean --force-pkgs-dirs --all --yes  && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> /etc/skel/.bashrc && \
    echo ". ${CONDA_DIR}/etc/profile.d/conda.sh && conda activate base" >> ~/.bashrc

ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]