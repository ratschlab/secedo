# splits BAM files by chromosome
step=10
base_dir="/cluster/work/grlab/projects/projects2019-supervario/simulated_data/varsim"
cov="cov01x"
n_cells=1000

mkdir -p "${base_dir}/${cov}/aligned_cells_split"
for idx in `seq 0 ${step} ${n_cells}`; do
  cmd1="echo hello"
  cmd2="echo hello"
  for i in `seq ${idx} $((idx+step-1))`; do
    suf=$(printf "%03d" $i)
    cmd1="$cmd1; ${base_dir}/split.sh healthy_${suf} ${base_dir}/${cov}"
    cmd2="$cmd2; ${base_dir}/split.sh tumor_${suf} ${base_dir}/${cov}"
  done
#  echo $cmd1
#  echo $cmd2
  bsub  -J "sp-${idx}" -W 1:00 -n 1 -R "rusage[mem=4000]" -R "span[hosts=1]" \
              -oo "${base_dir}/logs/sp-${idx}.lsf.log" "bash $cmd2; $cmd1"
done
