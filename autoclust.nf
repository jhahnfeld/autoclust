#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import java.nio.file.*

params.srv = 30
params.selector = 'eff'
params.id = params.srv
params.cov = params.srv
params.qcov = params.cov
params.scov = params.cov
params.minsize = 3
params.block = 100.KB

compressedFasta = params.input.endsWith(".gz")

println("Executor=" + params.executor)
println("Threads=" + params.threads)
println("SRV=" + params.srv)
params.mamba ? println("Using Mamba.") : println("Using Conda.")
println("\n")


workflow {
    fasta = Channel.fromPath(params.input)
    inflation = Channel.of(
            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
            3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0
          )

    preprocess_fasta(fasta)

    fasta_chunks = preprocess_fasta.out.faa.splitFasta(size: params.block, file: true, compress: true)
    fasta_count = preprocess_fasta.out.count.map{f -> f.text.toInteger()}

    makeblastdb(preprocess_fasta.out.faa)
    blast_input = fasta_chunks.combine(makeblastdb.out)
    blastp(blast_input, fasta_count)
    srv_transform(blastp.out)

    // save raw blast output
    compress(blastp.out.collectFile(name: 'blast.tsv',
                                    keepHeader: false,
                                    cache: 'lenient',
                                    sort: 'hash'))

    // prune singletons, pairs and nodes that are connected only to one node of a graph
    prune(srv_transform.out.collectFile(name: 'srv.tsv',
                                        keepHeader: false,
                                        cache: 'lenient',
                                        sort: 'hash'))

    // load to binary mcl matrix; make the matrix symmetric (mean of SRVs)
    mcxload(prune.out)
    // select k best edges of a node without introducing singletons
    mcxalter(mcxload.out)
    // write the pruned, symmetric and knn reduced data back to a abc TSV file
    mcxdump_abc(mcxalter.out)
    // separate unlinked subgraphs with similar attributes into separate files
    separate(mcxdump_abc.out)

    // load separated abc TSVs to binary mcl matrices
    // input: id_subgraph, abc TSV
    // ouput: id_subgraph, mcx, tab
    mcxload2(separate.out.flatten().map{[it.getBaseName().tokenize('.')[0].toInteger(), it]})

    // mcl clustering for all separated graphs and inflations
    // input: inflation, id_subgraph, mcx, tab
    // output: id_subgraph, inflation, mci.I
    mcl(inflation.combine(mcxload2.out))

    // look up the optimal inflation for the separated subgrpahs with clm info based on modularity
    // input: id_subgraph, []inflation, ...], [mci.I, ...], mcx, tab
    // output: id_subgraph, (best) inflation, (best) id_subgraph.mci.I, tab
    clm(mcl.out.groupTuple().join(mcxload2.out))

    // dump the clusters for the separated subgrpahs
    // input: id_subgraph, inflation, mci.I, tab
    // output: id_subgraph, inflation, cluster
    mcxdump_clusters(clm.out.cluster)

    // collect all clusters of the separated subgrpahs into one file
    // input: id_subgraph, inflation, clusters -> clusters
    // output: all clusters in one file
    autoclust = mcxdump_clusters.out.map{it[2]}.collectFile(name: 'auto.cluster.tsv',
                                                            keepHeader: false,
                                                            cache: 'lenient',
                                                            sort: 'hash')

    // group mci.I and tab files together by their subgraph id for all inflation values
    // input: id_subgraph, inflation, mci.I; id_subgraph, mcx, tab
    // output: inflation, mci.I, tab
    clusters_by_inflation = mcl.out.groupTuple().join(mcxload2.out).map{
        items -> {
        def combined = []
            if (items[1] instanceof List) {
                items[1].eachWithIndex{it, i -> {  //
                    combined.add([it,
                                  items[2][i],
                                  items[4]])}
                }
            } else {
                combined.add([items[1],
                              items[2],
                              items[4]])
            }
            return combined
        }
    }.flatten().collate(3)

    // dump clusters for all inflation values for each subgraph
    // input: inflation, mci.I, tab
    // output: inflation, clusters
    // mcxdump_for_cytoscape(clusters_by_inflation)

    // group clusters by inflation value
    //clusters_by_inflation = mcxdump_for_cytoscape.out.groupTuple().map{it[1]}.flatten().collectFile(
    //    keepHeader: false,
    //    cache: 'lenient',
    //    sort: 'hash'
    //).collect()

    // export the graph with cluster information for cytoscape
    // input: srv.tsv; autoclust_clusters; clusters_for_all_inflations
    // output: grraph_with_cluster_information
    //cytoscape(mcxdump_abc.out,
    //          autoclust,
    //          clusters_by_inflation)

    // export all clusters to fasta files
    // input: id_subgraph, inflation, clusters
    // output: id_subgraph.cluster_size.counter.fasta
    export_clusters_to_fasta(mcxdump_clusters.out)


    // add subgraph id and the number of sequences for each fasta file
    // input: id_subgraph.cluster_size.counter.fasta
    // output: id_subgraph.cluster_size.counter, sequence_count, id_subgraph.cluster_size.counter.fasta
    cluster_fasta = export_clusters_to_fasta.out.map{items -> {
                                                        def combined = []
                                                        if (items instanceof List) {
                                                            items.each{combined.add([it.getBaseName(),  // id_subgraph.cluster_size.counter
                                                                                     it.getBaseName().tokenize('.')[1].toInteger(),  // size
                                                                                     it])}  // path
                                                        } else {
                                                            combined.add([items.getBaseName(),  // id_subgraph.cluster_size.counter
                                                                          items.getBaseName().tokenize('.')[1].toInteger(),  // size
                                                                          items])  // path
                                                        }
                                                        return combined
                                                    }}.flatten().collate(3)

    // sort the fasta to the according muscle algorithm depending on their number of sequences
    aln = cluster_fasta.branch{
        ppp: it[1] <= 1000
        super5: it[1] > 1000
    }

    muscle_ppp(aln.ppp)
    muscle_super5(aln.super5)

    // build HMMs based on the multiple alignments
    // input: id_subgraph.cluster_size.counter, sequence_count, id_subgraph.cluster_size.counter.fasta, sprot.fasta
    // output: aln.id_subgraph.cluster_size.counter.afa
    hmm_build(muscle_ppp.out.concat(muscle_super5.out).collect(), fasta)
    compress_hmm(hmm_build.out.hmm)
}


process preprocess_fasta {
    cpus (params.threads >= 4 ? 4 : params.threads)
    memory { 4.GB * task.attempt }
    maxRetries 3

    input:
        path(faa)

    output:
        path('sorted.faa.gz'), emit: faa
        path('faa.count'), emit: count

    script:
    if ( compressedFasta == true )
        """
        pigz -dck $faa | grep -v '>' | sort -u | awk '{print ">"\$1"\\n"\$1}' | pigz -9 -p ${task.cpus} | \
        tee sorted.faa.gz | zgrep -c '>' > faa.count
        """
    else
        """
        grep -v '>' $faa | sort -u | awk '{print ">"\$1"\\n"\$1}' | pigz -9 -p ${task.cpus} | \
        tee sorted.faa.gz | zgrep -c '>' > faa.count
        """
}


process makeblastdb {
    memory { 4.GB * task.attempt }

    input:
        path(faa)

    output:
        path('blastdb/')

    script:
    """
    mkdir blastdb
    gzip -dck $faa | \
    makeblastdb -dbtype prot -hash_index -max_file_sz 4GB -title autoclust -out blastdb/db
    """
}


process blastp {
    memory { 8.GB * task.attempt }
    maxRetries 3
    array params.array

    input:
        tuple path(faa), path(db)
        val(count)

    output:
        path('blast.tsv')

    script:
    """
    pigz -dck $faa | \
    blastp -task blastp-short \
    -db $db/db \
    -out blast.tsv \
    -max_target_seqs $count \
    -matrix BLOSUM62 \
    -comp_based_stats 0 \
    -soft_masking false \
    -seg no \
    -evalue 10 \
    -outfmt '6 qseqid sseqid pident length bitscore evalue' \
    -num_threads ${task.cpus}
    """
}


process compress {
    publishDir "${params.output}", mode: 'copy', pattern: "*.gz"
    cpus (params.threads >= 8 ? 8: params.threads)
    memory { 1.GB * task.attempt }

    input:
        path('blast.tsv')

    output:
        path("*.gz")

    script:
    """
    pigz -f -p ${task.cpus} -9 blast.tsv
    """
}


process srv_transform {
    cpus (params.threads >= 2 ? 2 : params.threads)
    memory { 10.GB * task.attempt }
    maxRetries 3
    array params.array

    input:
        path('blast.tsv')

    output:
        path('srv.tsv')

    script:
    """
    blast_to_srv_matrix.py --input blast.tsv --mode blast --minsrv ${params.srv} --id ${params.id} \
    --query_cover ${params.qcov} --subject_cover ${params.scov} --threads ${task.cpus}
    """
}

process prune {
    input:
        path('srv.tsv')

    output:
        path('pruned.srv.tsv')

    script:
    """
    srv_pruning.py --input srv.tsv --low-memory --threads ${task.cpus}
    """
}


process mcxload {
    memory { 32.GB * task.attempt }

    input:
        path('srv.tsv')

    output:
        tuple path('srv.mcx'), path('srv.tab')

    script:
    """
    mcxload -ri add -tf 'scale(2),gt(${params.srv})' --write-binary -abc srv.tsv -write-tab srv.tab -o srv.mcx
    """
}
/*
    The square of the two SRVs of a protein pair might have interisting features in terms of similarity scaling.
    Use SRV * 10 as input. Squares are between 0 and 100.

    MINIMUM=\$((${params.srv}*${params.srv}/100))
    mcxload -ri mul -tf "scale(100),gt(\$MINIMUM)" --write-binary -abc srv.tsv -write-tab srv.tab -o srv.mcx
*/


process mcxalter {
    memory { 32.GB * task.attempt }

    input:
        tuple path('srv.unaltered.mcx'), path('srv.tab')

    output:
        tuple path('srv.mcx'), path('srv.tab')

    script:  // TODO add an additional step for k single digit selection
    """
    K=\$(top_k_edges.py --imx srv.unaltered.mcx -t ${task.cpus})
    echo \$K | tee k.txt

    mcx alter -imx srv.unaltered.mcx -tf "#knn(\$K)" --write-binary -o srv.mcx
    """
}


process mcxdump_abc {
    publishDir "${params.output}/srv", mode: 'copy', pattern: "pruned.srv.tsv"

    input:
        tuple path('srv.mcx'), path('srv.tab')

    output:
        path('pruned.srv.tsv')

    script:
    """
    mcxdump -imx srv.mcx -tab srv.tab -o pruned.srv.tsv
    """
}


process separate {
    memory { 16.GB * task.attempt }

    input:
        path('srv.tsv')

    output:
        path('*.subgraph.tsv')

    script:
    """
    separate_graph.py --abc srv.tsv --threads ${task.cpus}
    """
}


process mcxload2 {
    memory { 8.GB * task.attempt }

    input:
        tuple val(id), path('srv.tsv')

    output:
        tuple val(id), path('srv.mcx'), path('srv.tab')

    script:
    """
    mcxload --write-binary -abc srv.tsv -write-tab srv.tab -o srv.mcx
    """
}


process mcl {
    cpus { params.threads >= 4 * task.attempt * task.attempt ? 4 * task.attempt * task.attempt : params.threads }
    memory { 1.GB * task.attempt }
    array params.array

    input:
        tuple val(inflation), val(id), path('srv.mcx'), path('srv.tab')

    output:
        tuple val(id), val(inflation), path("*.clusters.mci.I*")

    script:
    """
    I=\$(echo $inflation | sed 's/\\.//g')
    echo \$I

    mcl srv.mcx -o ${id}.clusters.mci.I\$I -scheme 7 -I $inflation -te ${task.cpus}
    """
}


process clm {
    publishDir "${params.output}/clm/", mode: 'copy', pattern: "*.txt"
    memory { 4.GB * task.attempt }
    array params.array

    input:
        tuple val(id), val(inflations), path(clusters), path('srv.mcx'), path('srv.tab')

    output:
        tuple val(id), env(INFLATION), path("${id}.mci"), path('srv.tab'), emit: cluster
        path("*.txt"), emit: log

    script:
    """
    clm info -o ${id}.clminfo.txt 'srv.mcx' $clusters
    INFLATION=\$(mcl_inflation.py --input ${id}.clminfo.txt --id $id --metric ${params.selector})
    """
}


process mcxdump_clusters {
    array params.array

    input:
        tuple val(id), val(inflation), path(mci), path('srv.tab')

    output:
        tuple val(id), val(inflation), path("${id}.cluster.tsv")

    script:
    """
    mcxdump -icl $mci -tabr srv.tab -o ${id}.cluster.tsv
    """
}


process mcxdump_for_cytoscape {
    array params.array

    input:
        tuple val(inflation), path(mci), path('srv.tab')

    output:
        tuple val(inflation), path("cluster.${inflation}.tsv")

    script:
    """
    mcxdump -icl $mci -tabr srv.tab -o cluster.${inflation}.tsv
    """
}


process cytoscape {
    publishDir "${params.output}/srv", mode: 'copy', pattern: "srv.clustered.tsv"
    memory { 64.GB * task.attempt }

    input:
        path('srv.tsv')
        path(autoclust)
        path(clusters)

    output:
        path('srv.clustered.tsv')

    script:
    """
    export_to_cytoscape.py --abc srv.tsv --autoclust $autoclust --cluster $clusters
    """
}


process export_clusters_to_fasta {
    array params.array

    input:
        tuple val(id), val(inflation), path(clusters)

    output:
        path('*.faa')

    script:
    """
    split_clusters_for_alignment.py --input $clusters --id $id
    """
}


process muscle_ppp {
    publishDir "${params.output}/aln", mode: 'copy', pattern: "aln.*.afa"
    cpus (params.threads >= 8 ? 8 : params.threads)
    memory { 2.GB * task.attempt }
    maxRetries = 5
    array params.array

    input:
        tuple val(identifier), val(size), path(fasta)

    output:
        path("aln.${identifier}.afa")

    script:
    """
    muscle -amino -align $fasta -output aln.${identifier}.afa -threads ${task.cpus}
    """
}


process muscle_super5 {
    publishDir "${params.output}/aln", mode: 'copy', pattern: "aln.*.afa"
    cpus (params.threads >= 16 ? 16 : params.threads)
    memory { 4.GB * task.attempt }
    maxRetries = 3
    array params.array

    input:
        tuple val(identifier), val(size), path(fasta)

    output:
        path("aln.${identifier}.afa")

    script:
    """
    muscle -amino -super5 $fasta -output aln.${identifier}.afa -threads ${task.cpus}
    """
}


process hmm_build {
    publishDir "${params.output}", mode: 'copy', pattern: "clusters.*.tsv.gz"
    memory { 64.GB * task.attempt }

    input:
        path(aln)
        path(fasta)

    output:
        path("*.hmm"), emit: hmm
        path("clusters.tsv.gz"), emit: clusters

    script:
    """
    alignement_to_hmm.py --proteins $fasta --alignments ./
    """
}

process compress_hmm {
    publishDir "${params.output}", mode: 'copy', pattern: "*.hmm.gz"
    memory { 8.GB * task.attempt }

    input:
        path(hmm)

    output:
        path("*.hmm.gz")

    script:
    """
    pigz -f -p ${task.cpus} -9 $hmm
    """
}
