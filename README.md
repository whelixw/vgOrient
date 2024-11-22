STEP 1:
    Install conda: https://docs.anaconda.com/miniconda/install/#quick-command-line-install
STEP 2:
    conda install -y conda-build
    conda install -y bioconda-utils
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
STEP 3:
    conda build meta.yaml
STEP 4:
    conda create -n vgOrient-env vgoreint --use-local
    conda activate vgOrient-env
STEP 5:
    jaccard_dit_wrapper.py  datasets/suina/*.fasta --vg_output_dir VG_OUTPUT_DIR --output OUTPUT_NAME --orientation --min_jaccard_init -w 512 -m 256
