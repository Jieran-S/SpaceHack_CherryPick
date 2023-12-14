import os

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")

rule all:
    input:
        config['results_dir'] + "/" + config['dataset']


rule libd_dlpfc:
    output:
        dir=directory(config['results_dir'] + "/libd_dlpfc")
    conda:
        GIT_DIR + "/data/libd_dlpfc/libd_dlpfc.yml"
    shell:
        "Rscript {GIT_DIR}/data/libd_dlpfc/libd_dlpfc.r -o {output.dir}"

rule xenium_mouse_brain_SergioSalas:
    output:
        dir=directory(config['results_dir'] + "/xenium_mouse_brain_SergioSalas")
    conda:
        GIT_DIR + "/data/xenium-mouse-brain-SergioSalas/environment.yml"
    shell:
        "python {GIT_DIR}/data/xenium-mouse-brain-SergioSalas/prepare_data.py -o {output.dir}"

rule abc_atlas_wmb_thalamus:
    output:
        dir=directory(config['results_dir'] + "/abc_atlas_wmb_thalamus")
    conda:
        GIT_DIR + "/data/abc_atlas_wmb_thalamus/abc_atlas_wmb_thalamus.yml"
    shell:
        "python {GIT_DIR}/data/abc_atlas_wmb_thalamus/abc_atlas_wmb_thalamus.py -o {output.dir}"

rule SEA_AD_data:
    output:
        dir=directory(config['results_dir'] + "/SEA_AD_data")
    conda:
        GIT_DIR + "/data/SEA_AD_data/SEA_AD_data.yml"
    shell:
        "python {GIT_DIR}/data/SEA_AD_data/SEA_AD_data.py -o {output.dir}"


rule sotip_simulation:
    output:
        dir=directory(config['results_dir'] + "/sotip_simulation")
    conda:
        GIT_DIR + "/data/sotip_simulation/sotip.yml"
    shell:
        "python {GIT_DIR}/data/sotip_simulation/sotip.py -o {output.dir}"
