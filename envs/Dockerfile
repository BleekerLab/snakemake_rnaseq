FROM ubuntu:18.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/envs/env/bin:$PATH

# Add environment file where dependencies are listed
ADD environment.yaml /tmp/environment.yml


RUN apt-get update --fix-missing && \
    apt-get install -y wget  

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc 

RUN conda env create -f /tmp/environment.yml

RUN echo "conda activate counts" >> ~/.bashrc

CMD [ "/bin/bash" ]