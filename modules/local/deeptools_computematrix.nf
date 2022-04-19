process DEEPTOOLS_COMPUTEMATRIX {
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::deeptools=3.5.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.1--py_0' :
        'quay.io/biocontainers/deeptools:3.5.1--py_0' }"

    input:
    path bed // output bedfiles from MACS2_CALLPEAK
    path bigwig // output bigwig from DEEPTOOLS_BAMCOVERAGE
    path local_remote_files // paths for local or remote files (not generated in the pipline) to be staged
                                  // The paths for these files have to be explicitly specified in "assets/heatmap_blueprint.yaml"
    tuple val(name), //
          val(bed_name), // string of all bed files' names which is mix of files staged by path(bed) and path(local_remote_files)
          val(bigwig_name) // string of all bigwig files' names which is mix of files staged by path(bigwig) and path(local_remote_files)
    output:
    path "*.mat.gz"    , emit: matrix
    path "*.mat.tab"   , emit: table
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    computeMatrix \\
        $args \\
        --regionsFileName $bed_name \\
        --scoreFileName $bigwig_name \\
        --outFileName ${name}.computeMatrix.mat.gz \\
        --outFileNameMatrix ${name}.computeMatrix.vals.mat.tab \\
        --numberOfProcessors $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed -e "s/computeMatrix //g")
    END_VERSIONS
    """
}
