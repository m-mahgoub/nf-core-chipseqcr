process DEEPTOOLS_PLOTHEATMAP {
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    tuple val(name), //
          path (matrix), //
          val(bed_labels_inLine_with_quotes), // string of labels for all regions "beds" as specified in "assets/heatmap_blueprint.yaml"
          val(bigwig_labels_inLine_with_quotes) // string of labels for bigwig samples as specified in "assets/heatmap_blueprint.yaml"

    output:
    path "${name}.*"    , emit: plot
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def format = task.ext.format ?: 'png'
    def height = task.ext.height ?: 27
    def width = task.ext.width ?: 9
    """
    plotHeatmap \\
        $args \\
        --matrixFile $matrix \\
        --regionsLabel $bed_labels_inLine_with_quotes \\
        --samplesLabel $bigwig_labels_inLine_with_quotes \\
        --plotFileFormat $format \\
        --outFileName ${name}.${format} \\
        --heatmapHeight $height \\
        --heatmapWidth $width \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(plotHeatmap --version | sed -e "s/plotHeatmap //g")
    END_VERSIONS
    """
}
