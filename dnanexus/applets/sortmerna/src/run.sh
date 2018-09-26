#!/bin/bash

main() {
    #
    # prepare and validate environment i.e. binaries, libraries, PATH
    #
    OUT_DIR=$HOME/out
    IDX_DIR=$HOME/idx
    REFS_DIR=$HOME/in/refs
    READS_DIR=$HOME/in/reads
    DBDIR=kvdb
    TASK=4
    UPLOAD_ROOT=/out

    echo "mkdir $REFS_DIR"
    mkdir -p $REFS_DIR
    echo "mkdir $READS_DIR"
    mkdir -p $READS_DIR
    echo "mkdir $OUT_DIR"
    mkdir -p $OUT_DIR
    echo "mkdir $IDX_DIR"
    mkdir -p $IDX_DIR

    # set PATH
    chmod u+x $HOME/bin/indexdb
    chmod u+x $HOME/bin/sortmerna
    export PATH=$HOME/bin:$PATH

    echo "PATH: $PATH"

    # check patchelf
    echo "apt-cache policy patchelf"
    apt-cache policy patchelf
    echo "which patchelf"
    which patchelf # /usr/local/bin/patchelf

    # check rocksdb
    echo "apt-cache policy librocksdb-dev"
    apt-cache policy librocksdb-dev

    # check samtools
    echo "apt-cache policy samtools"
    apt-cache policy samtools

    # patch sortmerna to look for libstdc++.so.6 in the ORIGIN directory
    echo "Running: patchelf --set-rpath '$ORIGIN' $HOME/bin/sortmerna"
    patchelf --set-rpath '$ORIGIN' $HOME/bin/sortmerna

    echo "REFS:            ${REFS[@]}"
    echo "READS:           ${READS}"
    echo "READS_GZ:        ${READS_GZ}"
    echo "SAM:             ${SAM}"
    echo "FASTX:           ${FASTX}"
    echo "BLAST:           ${BLAST}"
    echo "advanced:        ${advanced}"
    echo "processing_task: ${TASK}"
    echo "help_smr:        ${help_smr}"

    # get executable name from the Job description to use in output (TODO: a better way?)
    tokens=($(dx describe ${DX_JOB_ID} | grep 'Executable name'))
    exe_name=${tokens[2]}
    upload_path="$UPLOAD_ROOT/$exe_name"
    echo "Applet/Executable name: $exe_name" # e.g. sortmerna-3.0-beta.run
    echo "Upload will be done into $upload_path"
    stat=$(dx mkdir -p $upload_path)
    echo "Upload path status: $stat"


    # print help strings if requested
    if [ "${help_smr}" == "true" ]; then
        indexdb -h || true
        echo
        echo
        sortmerna -h || true
    fi

    #
    # Download input files, process Options and run
    #
    #dx-download-all-inputs
    for REF in "${REFS[@]}"
    do
        echo "dx download $REFS_DIR/: \"$REF\""
        dx download -o $REFS_DIR/ "$REF"
    done

    reads_base_noext=""
    reads_ext=""

    if [[ ! -z "${READS}" ]]; then
        echo "download -o $READS_DIR/ \"${READS}\""
        dx download -o $READS_DIR/ "${READS}"
        reads_base_noext=${READS_name%.*} # strip extension
        reads_ext=${READS_name#*.} # strip name
        echo "reads_input basename: $reads_base_noext extension: $reads_ext"
    fi

    if [[ ! -z "${READS_GZ}" ]]; then
        echo "download -o $READS_DIR/ \"${READS_GZ}\""
        dx download -o $READS_DIR/ "${READS_GZ}"
        reads_base_noext=${READS_GZ_name%%.*} # strip extension
        tmp=${READS_GZ_name%.*} # strip 'gz'
        reads_ext=${tmp#*.}  # strip name
        echo "reads_input basename: $reads_base_noext extension: $reads_ext"
    fi

    echo "find $HOME -type f"
    find $HOME -type f

    # options
    opt_REF=""
    opt_READS=""
    opt_ALIGNED="--aligned $OUT_DIR/${reads_base_noext}_aligned"
    opt_OTHER="--other $OUT_DIR/${reads_base_noext}_other"
    opt_SAM=""
    opt_FASTX=""
    opt_LOG="--log"
    opt_BLAST=""
    opt_v_VERBOSE="-v"
    opt_d_KVDB="-d $DBDIR"
    opt_TASK="--task ${TASK}"
    #opt_a_NUM_THREADS="-a $(nproc)" # default anyway - ignore
    opts=""

    refopts=""
    num_refs=${#REFS[*]}
    for (( i=0; i<$(( num_refs )); ++i ))
    do
        bname=${REFS_name[$i]}
        refopt="$REFS_DIR/$bname,$IDX_DIR/${bname%.*}"
        if [[ ! -z "${refopts}" ]]; then
            refopts="${refopts}:${refopt}"
        else
            refopts=${refopt}
        fi
    done
    opt_REF="--ref ${refopts}"

    if [[ ! -z "${READS}" ]]; then
        opt_READS="--reads $READS_DIR/${READS_name}"
    fi

    if [[ ! -z "${READS_GZ}" ]]; then
    # TODO: update when SortMeRNA 2.2 is released (with reading gzip as option)
        opt_READS="--reads-gz $READS_DIR/${READS_GZ_name}"
    fi

    if [ "${SAM}" == "true" ]; then
        opt_SAM="--sam"
    fi

    if [ "${FASTX}" == "true" ]; then
        opt_FASTX="--fastx"
    fi

    # Build sortmerna index
    echo "Starting: indexdb $opt_REF $opt_v_VERBOSE"
    indexdb $opt_REF $opt_v_VERBOSE

    echo "Listing index: $IDX_DIR/ ..."
    ls -lrt $IDX_DIR

    # Run sortmerna
    if [[ ! -z "${BLAST}" ]]; then
        opts="$opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG --blast \"${BLAST}\" $advanced $opt_d_KVDB $opt_TASK"
        echo "Starting: sortmerna $opts"
        sortmerna $opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG --blast "${BLAST}" $advanced $opt_d_KVDB $opt_TASK
    else
        opts="$opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG $advanced $opt_d_KVDB $opt_TASK"
        echo "Starting: sortmerna $opts"
        sortmerna $opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG $advanced $opt_d_KVDB $opt_TASK
    fi

    #
    # print statistics
    #
    echo "[INFO] ls -lrt $OUT_DIR"
    ls -lrt $OUT_DIR

    # Log
    echo "[INFO] find $OUT_DIR -name '*.log' | xargs cat"
    find $OUT_DIR -name '*.log' | xargs cat

    # fasta line count
    echo "[INFO] find $OUT_DIR -name '*.fasta' | xargs wc -l"
    find $OUT_DIR -name '*.fasta' | xargs wc -l

    # fastq line count
    echo "[INFO] find $OUT_DIR -name '*.fastq' | xargs wc -l"
    find $OUT_DIR -name '*.fastq' | xargs wc -l

    # blast line count
    echo "[INFO] find $OUT_DIR -name '*.blast' | xargs wc -l"
    find $OUT_DIR -name '*.blast' | xargs wc -l

    # sam line count
    echo "[INFO] find $OUT_DIR -name '*.sam' | xargs wc -l"
    find $OUT_DIR -name '*.sqm' | xargs wc -l

    #"${reads_base_noext}_aligned.log"
    #"${reads_base_noext}_other.fasta"
    #"${reads_base_noext}_aligned.fasta"
    #"${reads_base_noext}_aligned.blast"

    #
    # Output results
    #
    # Aligned reads - required output
    mkdir -p $OUT_DIR/output_fastx_gz
    aligned_file="$OUT_DIR/output_fastx_gz/${reads_base_noext}_aligned.${reads_ext}.gz"
    echo "[INFO] Output aligned reads: $aligned_file"
    gzip -c $OUT_DIR/${reads_base_noext}_aligned.$reads_ext > $aligned_file
    file_id=$(dx upload $aligned_file --path "$upload_path/" --brief)
    dx-jobutil-add-output "output_fastx_gz" "$file_id"

    # Log file - required output
    mkdir -p $OUT_DIR/output_logfile
    log_file="$OUT_DIR/output_logfile/${reads_base_noext}_aligned.log"
    echo "[INFO] Output log file: $log_file"
    if [ -e "$OUT_DIR/${reads_base_noext}_aligned.log" ]; then
        mv $OUT_DIR/${reads_base_noext}_aligned.log $OUT_DIR/output_logfile/
    fi

    file_id=$(dx upload $log_file --path "$upload_path/" --brief)
    dx-jobutil-add-output "output_logfile" "$file_id"

    # Non-aligned reads
    mkdir -p $OUT_DIR/output_other_gz
    file_other_gz="$OUT_DIR/output_other_gz/${reads_base_noext}_other.${reads_ext}.gz"
    echo "[INFO] Output non-aligned reads: $file_other_gz"
    gzip -c $OUT_DIR/${reads_base_noext}_other.$reads_ext > $file_other_gz
    file_id=$(dx upload $file_other_gz --path "$upload_path/${reads_base_noext}_other.$reads_ext.gz" --brief)
    dx-jobutil-add-output "output_other_gz" "$file_id"

    # SAM (optional)
    if [ "${SAM}" == "true" ]; then
        echo "[INFO] Output BAM: $OUT_DIR/${reads_base_noext}_aligned.bam"

        refopts=""
        for (( i=0; i<$(( num_refs )); ++i ))
        do
            echo "[INFO] samtools faidx ${REFS_path[$i]}"
            samtools faidx ${REFS_path[$i]}
            refopts="${refopts} -t ${REFS_path[$i]}.fai"
        done
        
        echo "[INFO] samtools view ${refopts} -b $OUT_DIR/${reads_base_noext}_aligned.sam -o $OUT_DIR/${reads_base_noext}_aligned.bam"
        samtools view ${refopts} -b $OUT_DIR/${reads_base_noext}_aligned.sam -o $OUT_DIR/${reads_base_noext}_aligned.bam
        file_id=$(dx upload $OUT_DIR/${reads_base_noext}_aligned.bam --path "$upload_path/" --brief)
        dx-jobutil-add-output "output_sam" "$file_id"
    fi

    # BLAST (optional)
    if [[ ! -z "${BLAST}" ]]; then
        mkdir $OUT_DIR/output_blast_gz
        file_blast_gz="$OUT_DIR/output_blast_gz/${reads_base_noext}_aligned.blast.gz"
        echo "[INFO] Output BLAST: $file_blast_gz"
        gzip -c $OUT_DIR/${reads_base_noext}_aligned.blast > $file_blast_gz
        file_id=$(dx upload $file_blast_gz --path "$upload_path/${reads_base_noext}_aligned.blast.gz" --brief)
        dx-jobutil-add-output "output_blast_gz" "$file_id"
    fi

    echo "==== DONE ===="
}
