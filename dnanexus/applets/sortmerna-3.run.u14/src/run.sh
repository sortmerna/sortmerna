#!/bin/bash

main() {

    echo "[INFO] Running in $PWD"

    # set PATH
    export PATH=$HOME/bin:$PATH
    chmod u+x $HOME/bin/indexdb
    chmod u+x $HOME/bin/sortmerna
    ls -lrt $HOME/bin
    echo "[INFO] PATH: $PATH"

    # prepare directory structure
    IDX_DIR=$HOME/idx
    OUT_DIR=$HOME/out
    mkdir $IDX_DIR $OUT_DIR # otherwise ERROR: the directory /home/dnanexus/idx for writing index not found

    # download inputs
    dx-download-all-inputs

    echo "[INFO] Files in $HOME/"
    find . -type f

    for i in ${reference_input_path[@]}
    do
        echo "reference_input: $i"
    done
    echo "reads_input:     ${reads_input_path}"
    echo "reads_input_gz:  ${reads_input_gz_path}"
    echo "bam_output:      ${bam_output}"
    echo "fastx_output:    ${fastx_output}"
    echo "blast_output:    ${blast_output}"
    echo "advanced:        ${advanced}"
    echo "database_dir:    ${database_dir}"
    echo "processing_task: ${processing_task}"
    echo "output_dir:      ${output_dir}"
    echo "help_smr:        ${help_smr}"
    echo "dryidx:          ${dryidx}"
    echo "drysmr:          ${drysmr}"

    # get executable name from the Job description to use in output (TODO: is there a better way?)
    tokens=($(dx describe ${DX_JOB_ID} | grep 'Executable name'))
    exe_name=${tokens[2]}
    upload_path="${output_dir}/$exe_name"
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

    if [ ! -d "$database_dir" ]; then
        echo "Not Found KVDB directory: $database_dir"
    fi

    reads_file_basename=""
    ext=""

    if [[ ! -z "${reads_input}" ]]; then
        reads_file_basename=${reads_input_name%.*}
        ext=${reads_input_name#*.}
        echo "reads_input basename: $reads_file_basename extension: $ext"
    fi
    if [[ ! -z "${reads_input_gz}" ]]; then
        #dx download "${reads_input_gz}" -o ${reads_input_gz_name}
        gzip -d -c ${reads_input_gz_name} > ${reads_input_gz_name%.*}
        reads_input_name=${reads_input_gz_name%.*}
        reads_file_basename=${reads_input_name%.*}
        ext=${reads_input_name#*.}
        echo "reads_input extension: $ext"
    fi

    #
    # Set up options
    #
    opt_REF=""
    opt_READS=""
    opt_ALIGNED="--aligned $OUT_DIR/${reads_file_basename}_aligned"
    opt_OTHER="--other $OUT_DIR/${reads_file_basename}_other"
    opt_SAM=""
    opt_FASTX=""
    opt_LOG="--log"
    opt_BLAST=""
    opt_v_VERBOSE="-v"
    opt_d_KVDB="-d ${database_dir}"
    opt_TASK="--task ${processing_task}"
    #opt_a_NUM_THREADS="-a $(nproc)" # default anyway - ignore
    opts=""

    refs=""
    num_refs=${#reference_input[*]}
    for (( i=0; i<$(( num_refs )); ++i ))
    do
        ref=${reference_input_path[$i]}
        bname=${reference_input_name[$i]}
        refopt="$ref,$IDX_DIR/${bname%.*}"
        if [[ ! -z "${refs}" ]]; then
            refs="${refs}:${refopt}"
        else
            refs=${refopt}
        fi
    done
    opt_REF="--ref ${refs}"

    if [[ ! -z "${reads_input}" ]]; then
        opt_READS="--reads ${reads_input_path}"
    fi
    if [[ ! -z "${reads_input_gz}" ]]; then
    # TODO: update when SortMeRNA 2.2 is released (with reading gzip as option)
        opt_READS="--reads ${reads_input_path}"
    fi

      if [ "${bam_output}" == "true" ]; then
        opt_SAM="--sam"
    fi

    if [ "${fastx_output}" == "true" ]; then
        opt_FASTX="--fastx"
    fi

    # Build sortmerna index
    echo "Starting: indexdb $opt_REF $opt_v_VERBOSE"
    if [ "${dryidx}" == "false" ]; then
        indexdb $opt_REF $opt_v_VERBOSE
    fi
    echo "Listing index: $IDX_DIR/ ..."
    ls -lrt $IDX_DIR

    # Run sortmerna
    if [[ ! -z "${blast_output}" ]]; then
        opts="$opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG $opt_BLAST --blast \"${blast_output}\" $advanced $opt_d_KVDB $opt_TASK"
        echo "Starting: sortmerna $opts"
        if [ "${drysmr}" == "false" ]; then
            sortmerna $opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG --blast "${blast_output}" $advanced $opt_d_KVDB $opt_TASK
        fi
    else
        opts="$opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG $advanced $opt_d_KVDB $opt_TASK"
        echo "Starting: sortmerna $opts"
        if [ "${drysmr}" == "false" ]; then
            sortmerna $opt_REF $opt_READS $opt_ALIGNED $opt_OTHER $opt_SAM $opt_FASTX $opt_LOG $advanced $opt_d_KVDB $opt_TASK
        fi
    fi

    echo "Listing OUT idir: $OUT_DIR/ ..."
    ls -lrt $OUT_DIR

    #
    # Output results
    #
    # Aligned reads - required output
    mkdir -p $OUT_DIR/output_fastx_gz
    aligned_file="$OUT_DIR/output_fastx_gz/${reads_file_basename}_aligned.$ext.gz"
    echo "[INFO] Output aligned reads: $aligned_file"
    if [ "${drysmr}" == "false" ]; then
        gzip -c $OUT_DIR/${reads_file_basename}_aligned.$ext > $aligned_file
    else
        # prepare dummy output
        touch $aligned_file
        ls -l $aligned_file
    fi
    #file_id=$(dx upload $aligned_file --path $upload_path/${reads_file_basename}_aligned.$ext.gz --brief)
    file_id=$(dx upload $aligned_file --path "$upload_path/" --brief)
    dx-jobutil-add-output "output_fastx_gz" "$file_id"

    # Log file - required output
    mkdir -p $OUT_DIR/output_logfile
    if [ -e "$OUT_DIR/${reads_file_basename}_aligned.log" ]; then
        mv $OUT_DIR/${reads_file_basename}_aligned.log $OUT_DIR/output_logfile/
    fi
    log_file="$OUT_DIR/output_logfile/${reads_file_basename}_aligned.log"
    echo "[INFO] Output log file: $log_file"
    if [ "${drysmr}" == "true" ]; then
        # generate dummy output
        touch $log_file
        ls -l $log_file
    fi
    #file_id=$(dx upload ${reads_file_basename}_aligned.log --path "$upload_path/${reads_file_basename}_aligned.log" --brief)
    file_id=$(dx upload $log_file --path "$upload_path/" --brief)
    dx-jobutil-add-output "output_logfile" "$file_id"

    # Non-aligned reads
    mkdir -p $OUT_DIR/output_other_gz
    file_other_gz=$OUT_DIR/output_other_gz/${reads_file_basename}_other.$ext.gz
    echo "[INFO] Output non-aligned reads: $file_other_gz"
    if [ "${drysmr}" == "false" ]; then
        gzip -c $OUT_DIR/${reads_file_basename}_other.$ext > $file_other_gz
        file_id=$(dx upload $file_other_gz --path "$upload_path/${reads_file_basename}_other.$ext.gz" --brief)
        dx-jobutil-add-output "output_other_gz" "$file_id"
    fi

    # BAM (optional)
    if [ "${sam_output}" == "true" ]; then
        echo "[INFO] Output BAM: $OUT_DIR/${reads_file_basename}_aligned.bam"
        if [ "${drysmr}" == "false" ]; then
            samtools view -b ${reads_file_basename}_aligned.sam > ${reads_file_basename}_aligned.bam
            file_id=$(dx upload ${reads_file_basename}_aligned.bam --path "$upload_path/${reads_file_basename}_aligned.bam" --brief)
            dx-jobutil-add-output "output_bam" "$file_id"
        fi
    fi

    # BLAST (optional)
    if [[ ! -z "${blast_output}" ]]; then
        mkdir $OUT_DIR/output_blast_gz
        file_blast_gz=$OUT_DIR/output_blast_gz/${reads_file_basename}_aligned.blast.gz
        echo "[INFO] Output BLAST: $file_blast_gz"
        if [ "${drysmr}" == "false" ]; then
            gzip -c $OUT_DIR/${reads_file_basename}_aligned.blast > $file_blast_gz
            file_id=$(dx upload $file_blast_gz --path "$upload_path/${reads_file_basename}_aligned.blast.gz" --brief)
            dx-jobutil-add-output "output_blast_gz" "$file_id"
        fi
    fi

    #dx-upload-all-outputs

    echo "------------------ job_output.json ----------------------"
    cat job_output.json
    echo
    echo

    echo "------------------ job_input.json -----------------------"
    cat job_input.json
    echo
    echo "==== DONE ===="
}
