source activate munk
snakemake \
    --configfile $1 \
    --latency-wait 60 \
    -j 300 \
    --cluster-config cluster.yml \
    --cluster \
        "stdbuf -o0 -e0 sbatch --qos={cluster.qos} \
         --ntasks-per-node={cluster.cpu} \
         --mem={cluster.mem} \
         --time={cluster.time}" \
    rand_networks
mkdir -p slurm_logs
mv slurm*.out slurm_logs/
mv slurm*.err slurm_logs/
