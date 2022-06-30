process dehumanise {
    publishDir "$output_dir", mode: params.publish_mode, pattern: "no_human/*/*.fastq.gz"

    input:
    tuple val(sample), path(fq)

    output:
    tuple val(sample), path("no_human/${sample}/${sample}.bam"), emit: bam
    tuple val(sample), path("no_human/${sample}/${sample}*.fastq.gz"), emit: fq
    tuple val("${sample}.full"), path("${sample}.full.bam"), emit: full_bam

    script:
    def in_fq = (fq.size() == 2) ? "in=${fq[0]} in2=${fq[1]}" : "in=${fq[0]}";
    def maxmem = task.memory.toGiga()
    """
    set -o pipefail
    mkdir -p no_human/${sample}
    ln -s ${params.decon_ref}
    bbmap.sh -Xmx${maxmem}g t=$task.cpus ${in_fq} outu=unmapped.sam outm=mapped.sam idfilter=${params.decon_minid}
    cp unmapped.sam full.sam
    samtools view mapped.sam >> full.sam
    samtools view -f 4 mapped.sam >> unmapped.sam
    samtools view -f 8 mapped.sam >> unmapped.sam

    samtools collate -@ $task.cpus -O unmapped.sam | samtools view -buh > unmapped.bam

    if [[ "\$?" -eq 0 ]];
    then
        samtools fastq -@ task.cpus -0 ${sample}_other.fastq.gz -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz unmapped.bam
        if  [[ "\$?" -eq 0 ]];
        then

            if [[ -z "\$(gzip -dc ${sample}_R1.fastq.gz | head -n 1)" ]];
            then
                if [[ ! -z "\$(gzip -dc ${sample}_other.fastq.gz | head -n 1)" ]];
                then
                    mv -v ${sample}_other.fastq.gz no_human/${sample}/${sample}_R1.fastq.gz;
                fi;
            else
                    mv -v ${sample}_R1.fastq.gz no_human/${sample}/;
                    if [[ ! -z "\$(gzip -dc ${sample}_R2.fastq.gz | head -n 1)" ]];
                    then
                        mv -v ${sample}_R2.fastq.gz no_human/${sample}/;
                    fi;
            fi;

            mv -v unmapped.bam no_human/${sample}/${sample}.bam
            ls -l *.fastq.gz
            ls -l no_human/${sample}/
            rm -rf *.fastq.gz

            samtools view -buh full.sam > ${sample}.full.bam
        fi;
    fi;

    """
}
