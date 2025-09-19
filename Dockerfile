################################################
# Builds the cathy-mm image from github
#
# Run the following command first:
#  export DOCKER_BUILDKIT=1
#
# Then build with:
#  docker build -t cathy-mm:x.y .
#
# To run the image with display support:
# 1. On the host terminal: xhost +local:docker
# 2. Then the docker run line, for example:
#  docker run -v /home/username/sharedir:/sharedir -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY --volume="$HOME/.Xauthority:/root/.Xauthority:rw" --net=host cathy-mm cathy peek /sharedir/example
################################################

FROM mambaorg/micromamba:1.5.8-alpine3.20
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

USER root
RUN apk add --no-cache \
    openssh git \
    tzdata mesa-gl \
    glib libsm-dev libxrender libxext-dev

USER $MAMBA_USER

WORKDIR /opt/dock-cat/

ARG MAMBA_DOCKERFILE_ACTIVATE=1 # For pip and python commands to work

RUN git clone https://github.com/WrightGroupSRI/catheter_utils.git catheter_utils
RUN pip install -e catheter_utils

RUN git clone https://github.com/WrightGroupSRI/catheter_ukf.git catheter_ukf
RUN pip install -e catheter_ukf

RUN git clone https://github.com/WrightGroupSRI/dicom_utils.git dicom_utils
RUN pip install -e dicom_utils

RUN git clone https://github.com/WrightGroupSRI/dicom_art.git dicom_art
RUN pip install -e dicom_art

RUN git clone https://github.com/WrightGroupSRI/get_gt.git get_gt
RUN pip install -e get_gt

RUN git clone https://github.com/WrightGroupSRI/cathy.git cathy
RUN pip install -e cathy

COPY <<EOF /home/$MAMBA_USER/.jupyter/jupyter_notebook_config.py
c = get_config()
### If you want to auto-save .html and .py versions of your notebook:
# modified from: https://github.com/ipython/ipython/issues/8009
import os
from subprocess import check_call

def post_save(model, os_path, contents_manager):
    """post-save hook for converting notebooks to .py scripts"""
    if model['type'] != 'notebook':
        return # only do this for notebooks
    d, fname = os.path.split(os_path)
    check_call(['jupyter', 'nbconvert', '--to', 'script', fname], cwd=d)
    check_call(['jupyter', 'nbconvert', '--to', 'html', fname], cwd=d)

c.FileContentsManager.post_save_hook = post_save
EOF

USER root
RUN chown -R $MAMBA_USER:$MAMBA_USER /home/$MAMBA_USER/.jupyter
RUN chmod a+w /home/$MAMBA_USER/.jupyter
USER $MAMBA_USER

