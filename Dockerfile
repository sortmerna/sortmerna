
FROM ubuntu:jammy

SHELL ["/bin/bash", "-c"]
RUN apt update && apt upgrade -y && \
    apt install -y build-essential git wget file

RUN wget -q -P /tmp https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    /bin/bash /tmp/Miniforge3-Linux-x86_64.sh -b -p /opt/miniforge3 && \
    rm /tmp/Miniforge3-Linux-x86_64.sh && \
    ln -s /opt/miniforge3/bin/conda /usr/local/bin && \
    ln -s /opt/miniforge3/bin/mamba /usr/local/bin && \
    mamba init && \
    mamba install -y pyyaml jinja2 requests ninja cmake=3.29.6
    # scikit-bio

ENV PATH=/opt/miniforge3/bin:$PATH
RUN git clone https://github.com/biocore/sortmerna.git && \
    cd /sortmerna && \
    python setup.py -n all




