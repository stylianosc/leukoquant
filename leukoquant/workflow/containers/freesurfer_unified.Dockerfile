FROM continuumio/miniconda3:25.3.1-1

# Build context must be the workflow/ directory:
#   docker buildx build --platform linux/amd64 \
#     -f containers/freesurfer_unified.Dockerfile -t ... .

WORKDIR /

# ── system utils ──────────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    bc libgomp1 perl tar tcsh wget vim-common \
    mesa-utils libxext6 libsm6 libxrender1 libxmu6 \
    libgl1 libglx-mesa0 \
    curl unzip file dc \
    pulseaudio libquadmath0 libgtk2.0-0 \
    gcc time && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# ── FreeSurfer 7.4.1 ──────────────────────────────────────────────────────────
ADD https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/7.4.1/freesurfer-linux-centos7_x86_64-7.4.1.tar.gz ./fs.tar.gz
RUN tar --no-same-owner -xzvf fs.tar.gz && \
    mv freesurfer /usr/local && \
    rm -f fs.tar.gz

ENV OS=Linux
ENV FREESURFER_HOME=/usr/local/freesurfer
ENV FREESURFER=/usr/local/freesurfer
ENV SUBJECTS_DIR=/usr/local/freesurfer/subjects
ENV LOCAL_DIR=/usr/local/freesurfer/local
ENV FSFAST_HOME=/usr/local/freesurfer/fsfast
ENV FMRI_ANALYSIS_DIR=/usr/local/freesurfer/fsfast
ENV FUNCTIONALS_DIR=/usr/local/freesurfer/sessions
ENV FS_OVERRIDE=0
ENV FIX_VERTEX_AREA=""
ENV FSF_OUTPUT_FORMAT=nii.gz
ENV PYTHONUNBUFFERED=1
ENV MINC_BIN_DIR=/usr/local/freesurfer/mni/bin
ENV MINC_LIB_DIR=/usr/local/freesurfer/mni/lib
ENV MNI_DIR=/usr/local/freesurfer/mni
ENV MNI_DATAPATH=/usr/local/freesurfer/mni/data
ENV MNI_PERL5LIB=/usr/local/freesurfer/mni/share/perl5
ENV PERL5LIB=/usr/local/freesurfer/mni/share/perl5
ENV PATH=$PATH:/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin

# ── MATLAB Runtime ────────────────────────────────────────────────────────────
RUN fs_install_mcr R2019b && \
    rm -rf /tmp/MCR* /tmp/matlab* 2>/dev/null || true

# ── FSL ───────────────────────────────────────────────────────────────────────
ENV FSLDIR="/usr/local/fsl"
ENV DEBIAN_FRONTEND="noninteractive"
ENV LANG="en_GB.UTF-8"

ADD https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/fslinstaller.py ./fslinstaller.py
RUN python ./fslinstaller.py -d /usr/local/fsl/ && \
    rm -f /fslinstaller.py
ENV PATH="/usr/local/fsl/bin:$PATH"

# ── ANTs 2.6.5 ────────────────────────────────────────────────────────────────
ADD https://github.com/ANTsX/ANTs/releases/download/v2.6.5/ants-2.6.5-centos7-X64-gcc.zip /ants.zip
RUN mkdir -p /ants && \
    unzip /ants.zip -d /ants && \
    rm -f /ants.zip
ENV ANTSPATH=/ants/ants-2.6.5/bin
ENV PATH=$ANTSPATH:$PATH

# ── NiftyReg ─────────────────────────────────────────────────────────────────
ADD https://github.com/KCL-BMEIS/niftyreg/releases/download/v2.0.0/NiftyReg-Ubuntu-v2.0.0.zip ./NiftyReg.zip
RUN unzip NiftyReg.zip -d /usr/local/niftyreg && \
    rm -f NiftyReg.zip
ENV PATH="/usr/local/niftyreg/NiftyReg-Ubuntu-v2.0.0:${PATH}"

# ── MRtrix3 ───────────────────────────────────────────────────────────────────
RUN conda install -c conda-forge -c MRtrix3 mrtrix3 libstdcxx-ng && \
    conda clean --all --force-pkgs-dirs -y

# ── Python packages (unified env baked into base) ─────────────────────────────
COPY envs/freesurfer_unified_env.yaml /tmp/freesurfer_unified_env.yaml
RUN conda env update --name base --file /tmp/freesurfer_unified_env.yaml && \
    conda clean --all --force-pkgs-dirs -y && \
    pip cache purge && \
    rm -f /tmp/freesurfer_unified_env.yaml && \
    find /opt/conda -type f -name '*.pyc' -delete && \
    find /opt/conda -type d -name '__pycache__' -exec rm -rf {} + 2>/dev/null || true

# ── Final sweep ───────────────────────────────────────────────────────────────
RUN rm -rf /tmp/* /var/tmp/* /root/.cache 2>/dev/null || true

ENTRYPOINT [ "sh", "-c", ". /usr/local/fsl/etc/fslconf/fsl.sh && /bin/bash" ]
