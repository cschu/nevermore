process run_samestr_convert {
    
    input:
		path(mp_sam)
        path(mp_profile)
		path(marker_db)

    output:
        path("sstr_convert/*/*.npz"), emit: sstr_npy

    script:
    """
    samestr --verbosity DEBUG \
    convert \
        --input-files ${mp_sam} \
        --min-vcov 1 \
        --min-aln-qual 0 \
        --mp-profiles-extension .txt \
        --marker-dir ${marker_db} \
        --output-dir sstr_convert/ \
        --nprocs ${task.cpus}
    """
}

process run_samestr_merge {
    
    input:
        tuple val(species), path(sstr_npy)

    output:
        tuple \
            val(species), \
            path("sstr_merge/${species}.npz"), \
            path("sstr_merge/${species}.names.txt"), \
        emit: sstr_npy

    script:
    """
    samestr --verbosity DEBUG \
    merge \
        --input-files ${sstr_npy} \
        --output-dir sstr_merge/ \
        --species ${species} \
        --nprocs ${task.cpus}
    """
}

process run_samestr_filter {
    
    input:
        tuple val(species), path(sstr_npy), path(sstr_names)
		path(marker_db)

    output:
        tuple \
            val(species), \
            path("sstr_filter/${species}.npz"), \
            path("sstr_filter/${species}.names.txt"), \
        emit: sstr_npy, optional: true

    script:
    // #    --global-pos-min-n-vcov 10 \
    // #    --sample-pos-min-n-vcov 2 \
    """

    samestr --verbosity DEBUG \
    filter \
        --input-files ${sstr_npy} \
        --input-names ${sstr_names} \
        --output-dir sstr_filter/ \
        --marker-dir ${marker_db} \
        --marker-trunc-len 50 \
        --global-pos-min-f-vcov 0.25 \
        --sample-pos-min-sd-vcov 3 \
        --samples-min-n-hcov 5000 \
        --sample-var-min-n-vcov 2 \
        --sample-var-min-f-vcov 0.025 \
        --species-min-samples 1 \
        --nprocs ${task.cpus}
    """
}

process run_samestr_stats {
    
    input:
        tuple val(species), path(sstr_npy), path(sstr_names)

    output:
        path "sstr_stats/${species}.aln_stats.txt", emit: sstr_stats

    script:
    """
    samestr --verbosity DEBUG \
    stats \
    --input-files ${sstr_npy} \
    --input-names ${sstr_names} \
    --nprocs ${task.cpus} \
    --output-dir sstr_stats/
    """
}

process run_samestr_compare {
    
    input:
        tuple val(species), path(sstr_npy), path(sstr_names)

    output:
        tuple \
            path("sstr_compare/${species}.closest.txt"), \
            path("sstr_compare/${species}.fraction.txt"), \
            path("sstr_compare/${species}.overlap.txt"), \
        emit: sstr_compare

    script:
    """
    samestr --verbosity DEBUG \
    compare \
        --input-files ${sstr_npy} \
        --input-names ${sstr_names} \
        --output-dir sstr_compare/ \
        --nprocs ${task.cpus}
    """
}

process run_samestr_summarize {
    
    input:
        path(sstr_data)
        path(mp_profiles)

    output:
        tuple \
            path("sstr_summarize/mp_counts.tsv"), \
            path("sstr_summarize/mp_species.tsv"), \
            path("sstr_summarize/mp_taxonomy.tsv"), \
            path("sstr_summarize/sstr_cooccurrences.tsv"), \
            path("sstr_summarize/sstr_strain_events.tsv"), \
        emit: sstr_summarize

    script:
    """
    mkdir profiles/
    find . -maxdepth 1 -name '*.mp4.txt' -exec mv -v {} profiles/ \\;

    samestr --verbosity DEBUG \
    summarize \
        --input-dir ./ \
        --mp-profiles-dir ./profiles/ \
        --mp-profiles-extension .txt \
        --output-dir sstr_summarize/
    """
}