run_iqtree() {
  local output_dir="$1"
  local fasta_file="$2"
  local model="$3"
  local n_cores="$4"
  local fast_str="$5"
  local pers_val="$6"
  local log_file="$7"
  local seed_local="$8"

  local output_pre="$output_dir/run_$seed_local"

  if [[ "$fast_str" == "true" ]]; then
    iqtree2 -s "$fasta_file" -m "$model" -nt "$n_cores" -seed "$seed_local" -pre "$output_pre" -pers "$pers_val" -safe -quiet -fast
  else
    iqtree2 -s "$fasta_file" -m "$model" -nt "$n_cores" -seed "$seed_local" -pre "$output_pre" -pers "$pers_val" -safe -quiet 
  fi

  echo "Finished iteration "$seed_local" at: $(date)" >> "$log_file"

  # save disk space by deleting files not needed for CellREST
  rm -f "$output_pre.parstree" "$output_pre.ckp.gz" "$output_pre.mldist" "$output_pre.uniqueseq.phy" "$output_pre.varsites.phy" "$output_pre.iqtree" "$output_pre.bionj"
}
