#source ~/miniconda3/bin/activate lukonda

if [ "$HOSTNAME" == "n104" ] || [ "$HOSTNAME" == "n105" ]; then
    echo 'Running on numbers cluster'
    snakemake -p \
        --cluster-config "cluster.json" \
        --drmaa ' --mem-per-cpu={cluster.mem} {cluster.flags} --nodes=1 --ntasks={cluster.threads} --nodes=1' \
        --jobs $1 \
        --latency-wait 60 \
        --max-jobs-per-second 5 \
        --restart-times 1 \
        --rerun-incomplete \
        --keep-going \
	-s $2 \
	${@:3}
    mkdir -p logs/cluster_logs
    mv slurm*.out logs/cluster_logs

elif [ "$HOSTNAME" == "login-apollo.hpc.bcgsc.ca" ]; then
    echo "Running on APOLLO cluster"
    mkdir -p logs/sge_logs_error logs/sge_logs_output
    snakemake -p \
        --cluster-config "config/cluster.json" \
        --cluster "qsub -q arc.q -P arc.prj -V -N msig_timing -e logs/sge_logs_error -o logs/sge_logs_output -l h_vmem={cluster.mem}M -pe smp {cluster.threads}" \
        --jobs 500 \
        --latency-wait 60 \
        --max-jobs-per-second 1 \
        --restart-times 5 \
        --rerun-incomplete

else
    echo "Running in non-cluster mode"

    snakemake -p \
        --cores $1 \
        --latency-wait 60 \
        --max-jobs-per-second 1 \
        --restart-times 5 \
        --rerun-incomplete \
        -s $2
fi
