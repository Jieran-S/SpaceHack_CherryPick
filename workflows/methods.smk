import os
from shared.functions import get_sample_dirs, get_ncluster

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")
TECHNOLOGY = config["technology"]
SEED = "42"

#####

def create_input_config(dataset):
    input_files = []
    sample_dirs = get_sample_dirs(config["data_dir"])
    for sample_dir in sample_dirs:
        for config_file_name in config["config_files"][dataset].keys():
            input_files.append(sample_dir + "/" + dataset + "/" + config_file_name + "/domains.tsv")
    return input_files

def create_input(dataset):
    input_files = []
    sample_dirs = get_sample_dirs(config["data_dir"])
    for sample_dir in sample_dirs:
        input_files.append(sample_dir + "/" + dataset + "/domains.tsv")
    return input_files

# without config
def create_BayesSpace_input(wildcards):
    return create_input("BayesSpace")

def create_scMEB_input(wildcards):
    return create_input("scMEB")

def create_SCAN_IT_input(wildcards):
    return create_input("SCAN_IT")

# with config
def create_spaGCN_input(wildcards):
    return create_input_config("spaGCN")

def create_GraphST_input(wildcards):
    return create_input_config("GraphST")

def create_BANKSY_input(wildcards):
    return create_input_config("BANKSY")

def create_meringue_input(wildcards):
    return create_input_config("meringue")

rule all:
    input: 
        create_BayesSpace_input,
        create_scMEB_input,
        #create_SCAN_IT_input, unklar
        create_spaGCN_input,
        #create_GraphST_input, unklar
        #create_BANKSY_input, unklar
        create_meringue_input

rule spaGCN:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        image = config["data_dir"] + "/{sample}/H_E.tiff",
    output:
        dir = directory(config["data_dir"] + "/{sample}/spaGCN/{config_file_name}"),
        file = config["data_dir"] + "/{sample}/spaGCN/{config_file_name}/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED,
        configfile = lambda wildcards: config["config_files"]["spaGCN"][wildcards.config_file_name]
    conda:
        GIT_DIR + "/method/spaGCN/spaGCN.yml"
    shell:
        """
        {GIT_DIR}/method/spaGCN/spaGCN.py \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            --image {input.image} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed} \
            --config {GIT_DIR}/method/spaGCN/{params.configfile}
        """


rule BayesSpace:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        dim_red = config["data_dir"] + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
    output:
        dir = directory(config["data_dir"] + "/{sample}/BayesSpace"),
        file = config["data_dir"] + "/{sample}/BayesSpace/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED
    conda:
        GIT_DIR + "/method/BayesSpace/BayesSpace.yml"
    shell:
        """
        {GIT_DIR}/method/BayesSpace/BayesSpace.r \
            -c {input.coordinates} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            --dim_red {input.dim_red} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed}
        """

rule scMEB_requirements:
    output:
        temp("scMEB_requirements.info")
    conda:
        GIT_DIR + "/method/SC.MEB/SC.MEB.yml"
    shell:
        """
        conda run R -e \"if(!require(SC.MEB)) install.packages('SC.MEB',repos = 'https://cran.r-project.org/')\"
        touch scMEB_requirements.info
        """

rule scMEB:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        dim_red = config["data_dir"] + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
        neighbors = config["data_dir"] + "/{sample}/delaunay_triangulation.mtx",
        requirements = "scMEB_requirements.info"
    output:
        dir = directory(config["data_dir"] + "/{sample}/scMEB"),
        domains = config["data_dir"] + "/{sample}/scMEB/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED
    conda:
        GIT_DIR + "/method/SC.MEB/SC.MEB.yml"
    shell:
        """
        /usr/bin/env Rscript {GIT_DIR}/method/SC.MEB/SC.MEB.r \
            -c {input.coordinates} \
            -f {input.features} \
            -o {input.observations} \
            -n {input.neighbors} \
            -d {output.dir} \
            --dim_red {input.dim_red} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed}
        """


rule GraphST:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
    output:
        dir = directory(config["data_dir"] + "/{sample}/GraphST/{config_file_name}"),
        file = config["data_dir"] + "/{sample}/GraphST/{config_file_name}/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED,
        configfile = lambda wildcards: config["config_files"]["GraphST"][wildcards.config_file_name]
    conda:
        GIT_DIR + "/method/GraphST/GraphST.yml"
    shell:
        """
        python {GIT_DIR}/method/GraphST/method_GraphST.py \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed} \
            --config {GIT_DIR}/method/spaGCN/{params.configfile}
        """


rule BANKSY_requirements:
    output:
        temp("BANKSY_requirements.info")
    conda:
        GIT_DIR + "/method/BANKSY/banksy.yml"
    shell:
        """
        conda run Rscript -e \"remotes::install_github('prabhakarlab/Banksy@v0.1.5', dependencies = TRUE)\"
        touch BANKSY_requirements.info
        """

rule BANKSY:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        requirements = "BANKSY_requirements.info"
    output:
        dir = directory(config["data_dir"] + "/{sample}/BANKSY/{config_file_name}"),
        file = config["data_dir"] + "/{sample}/BANKSY/{config_file_name}/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED,
        configfile = lambda wildcards: config["config_files"]["BANKSY"][wildcards.config_file_name]
    conda:
        GIT_DIR + "/method/BANKSY/banksy.yml"
    shell:
        """
        Rscript {GIT_DIR}/method/BANKSY/banksy.r \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed} \
            --config {GIT_DIR}/method/BANKSY/{params.configfile}
        """

rule meringue_requirements:
    output:
        temp("meringue_requirements.info")
    conda:
        GIT_DIR + "/method/meringue/meringue.yml"
    shell:
        """
        conda run Rscript -e \"remotes::install_github('JEFworks-Lab/MERINGUE', ref = 'ca9e2ccabd95680d9ca0b323a8a507c038f2ea13')\"
        touch meringue_requirements.info
        """

rule meringue:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        neighbors = config["data_dir"] + "/{sample}/delaunay_triangulation.mtx",
        dim_red = config["data_dir"] + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
        requirements = "meringue_requirements.info",
    output:
        dir = directory(config["data_dir"] + "/{sample}/meringue/{config_file_name}"),
        file = config["data_dir"] + "/{sample}/meringue/{config_file_name}/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED,
        configfile = lambda wildcards: config["config_files"]["meringue"][wildcards.config_file_name]
    conda:
        GIT_DIR + "/method/meringue/meringue.yml"
    shell:
        """
        Rscript {GIT_DIR}/method/meringue/meringue.r \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -n {input.neighbors} \
            -d {output.dir} \
            --n_clusters {params.n_clusters} \
            --dim_red {input.dim_red} \
            --technology {params.technology} \
            --seed {params.seed} \
            --config {GIT_DIR}/method/meringue/{params.configfile}
        """

rule SCAN_IT:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        neighbors = config["data_dir"] + "/{sample}/delaunay_triangulation.mtx",
        dim_red = config["data_dir"] + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
    output:
        dir = directory(config["data_dir"] + "/{sample}/SCAN_IT/"),
        file = config["data_dir"] + "/{sample}/SCAN_IT/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED,
        configfile = GIT_DIR + "/method/SCAN-IT/config.json"
    conda:
        GIT_DIR + "/method/SCAN-IT/scanit.yml"
    shell:
        """
        python {GIT_DIR}/method/SCAN-IT/method_scanit.py \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -n {input.neighbors} \
            -d {output.dir} \
            --n_clusters {params.n_clusters} \
            --dim_red {input.dim_red} \
            --technology {params.technology} \
            --seed {params.seed} \
            --config {GIT_DIR}/method/meringue/{params.configfile}
        """