#!/bin/bash
# stripes/stripenn/scripts/math.sh
# Verify 250831 merged label swap: .hic fixed, .mcool still swapped

ctrl_M1=580458938
ctrl_M2=554893997
ctrl_M3=595748617
mut_M1=600982570
mut_M2=602015250
mut_M3=611092030

ctrl_sum=$(( ctrl_M1 + ctrl_M2 + ctrl_M3 ))
mut_sum=$(( mut_M1 + mut_M2 + mut_M3 ))

hic_ctrl=1731101552
hic_mut=1814089850
mcool_ctrl=1814089850
mcool_mut=1731101552

printf "ctrl reps sum: %d\n" "${ctrl_sum}"
printf "mut  reps sum: %d\n" "${mut_sum}"
printf "\n"

if [ "${hic_ctrl}" -eq "${ctrl_sum}" ]; then
  printf ".hic  ctrl_merged (%d) == ctrl reps sum (%d) => YES\n" "${hic_ctrl}" "${ctrl_sum}"
else
  printf ".hic  ctrl_merged (%d) != ctrl reps sum (%d) => NO\n" "${hic_ctrl}" "${ctrl_sum}"
fi

if [ "${hic_mut}" -eq "${mut_sum}" ]; then
  printf ".hic  mut_merged  (%d) == mut reps sum  (%d) => YES\n" "${hic_mut}" "${mut_sum}"
else
  printf ".hic  mut_merged  (%d) != mut reps sum  (%d) => NO\n" "${hic_mut}" "${mut_sum}"
fi

printf "\n"

if [ "${mcool_ctrl}" -eq "${mut_sum}" ]; then
  printf ".mcool ctrl_merged (%d) == MUT reps sum (%d) => SWAPPED\n" "${mcool_ctrl}" "${mut_sum}"
else
  printf ".mcool ctrl_merged (%d) != mut reps sum (%d) => ok\n" "${mcool_ctrl}" "${mut_sum}"
fi

if [ "${mcool_mut}" -eq "${ctrl_sum}" ]; then
  printf ".mcool mut_merged  (%d) == CTRL reps sum (%d) => SWAPPED\n" "${mcool_mut}" "${ctrl_sum}"
else
  printf ".mcool mut_merged  (%d) != ctrl reps sum (%d) => ok\n" "${mcool_mut}" "${ctrl_sum}"
fi

printf "\nConclusion: .hic fixed, .mcool needs reconversion (Stage 0)\n"
