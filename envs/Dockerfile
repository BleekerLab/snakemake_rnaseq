FROM ubuntu:18.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH=/opt/conda/bin:$PATH
ENV PATH=/opt/conda/envs/env/bin:$PATH

###############
# General setup
###############

# Create a non-root user
ENV username=al-khawarizmi
ENV uid=1000
ENV gid=100
ARG USERNAME=marc
ENV UID=$uid
ENV GID=$gid
ENV HOME=/home/$USER

# # Change Shell as conda is not compatible with default /bin/sh
SHELL [ "/bin/bash", "--login", "-c" ]


RUN adduser --home $HOME --disabled-password --gecos "Non-root-user" --uid $UID --gid $GID $USERNAME
    

# Add Miniconda3 script
ADD Miniconda3-py37_4.8.3-Linux-x86_64.sh /tmp/Miniconda3-py37_4.8.3-Linux-x86_64.sh  
# Add environment file where dependencies are listed
ADD environment.yaml /tmp/environment.yaml

# Add entrypoint script to be executed so that conda env is activated
ADD entrypoint.sh /usr/local/bin/


RUN chown $UID:$GID /usr/local/bin/entrypoint.sh 
RUN chmod u+x /usr/local/bin/entrypoint.sh

# # wget is required to download miniconda installer script
RUN apt-get update --fix-missing &&\
    apt-get install -y wget  

# # create a project directory inside user home
ENV PROJECT_DIR $HOME/rnaseq
RUN mkdir $PROJECT_DIR
WORKDIR $PROJECT_DIR


USER $USERNAME
####################
## install miniconda
####################

# # download and install Miniconda
ENV CONDA_DIR $HOME/miniconda3
RUN chmod u+x /tmp/Miniconda3-py37_4.8.3-Linux-x86_64.sh
RUN bash /tmp/Miniconda3-py37_4.8.3-Linux-x86_64.sh -b -p $CONDA_DIR
RUN rm /tmp/Miniconda3-py37_4.8.3-Linux-x86_64.sh

# # make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# # make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

# # make conda activate command available from /bin/bash --interative shells
RUN conda init bash


###########################
## Create conda env
###########################

RUN conda env create --name rnaseq --file /tmp/environment.yaml --force 
RUN conda clean --all --yes

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]
