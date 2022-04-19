process HEATMAP_BLUEPRINT {
    label 'process_very_low'

    conda (params.enable_conda ? 'bioconda::pyyaml=6.0 bioconda::pandas=1.4.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'mdivr/conda-nf-cutnrun:v0.1' :
        'mdivr/conda-nf-cutnrun:v0.1' }"

    input:
        path(heatmap_blueprint) // user defined file for plotting strategy in YAML format in “assets/heatmap_blueprint.yaml”
    output:
    path 'names_labels_for_deeptools.txt' , emit: str_for_deeptool
    path 'local_remote_paths_for_deeptools.txt' , emit: path_for_deeptool
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    heatmap_blueprint.py \\
    $heatmap_blueprint \\
    names_labels_for_deeptools.txt \\
    local_remote_paths_for_deeptools.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
