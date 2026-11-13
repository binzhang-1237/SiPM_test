#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
input_file="${1:-$script_dir/data/Run94_list.csv}"
output_dir="${2:-$script_dir/plots}"

# Write build output, peak-fit results, and ROOT plotting messages to one log file.
mkdir -p "$script_dir/log"
log_file="$script_dir/log/$(basename -- "${input_file%.*}").log"
echo "Log: $log_file"
exec >"$log_file" 2>&1

echo "Started: $(date --iso-8601=seconds)"
echo "Input: $input_file"
echo "Output: $output_dir"
bash "$script_dir/build.sh"

exec "$script_dir/build/analyze_sipm" "$input_file" "$output_dir"
