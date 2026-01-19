# This script is meant to convert all conditions specificied in the well_specification.csv file to a matched comparison file
## We will just look for a reference condition first (DMSO, WT, etc), and put it in the first column, and then look for all other conditions to compare against it
## The output comparison file will have the format:
## reference condition <tab> comparison condition1 <tab> comparison condition2 ...
#!/usr/bin/env bash

well_specification_file="E2904_well_specification.csv"
comparison_file="mt20260112_comparisons.txt"

reference_condition=$(grep -E "DMSO|WT" "$well_specification_file" | cut -d',' -f2 | head -1)
if [ -z "$reference_condition" ]; then
  reference_condition=$(cut -d',' -f2 "$well_specification_file" | head -1)
  echo "No DMSO/WT found; using first condition as reference: $reference_condition"
else
  echo "Reference condition identified: $reference_condition"
fi

# Remove any stray carriage returns from the reference condition
reference_condition=${reference_condition//$'\r'/}

mapfile -t all_conditions < <(awk -F',' '{print $2}' "$well_specification_file" | sed 's/\r$//' | awk '!seen[$0]++')

others=()
for c in "${all_conditions[@]}"; do
  if [ "$c" != "$reference_condition" ]; then
    others+=("$c")
  fi
done

if [ ${#others[@]} -eq 0 ]; then
  printf "%s\n" "$reference_condition" > "$comparison_file"
else
  # Join with tab using printf and remove any carriage returns to ensure a single-row output
  others_line=$(printf "%s\t" "${others[@]}")
  others_line=${others_line%$'\t'}
  others_line=${others_line//$'\r'/}
  printf "%s\t%s\n" "$reference_condition" "$others_line" > "$comparison_file"
fi

echo "Comparison file generated: $comparison_file"