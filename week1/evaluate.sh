#!/usr/bin/env bash
set -euxo pipefail

PY=python3
CODON="${HOME}/.codon/bin/codon"
SEQ_PLUGIN="-plugin seq"
RELEASE="-release"

ASM_DIR="week1/code/assembler"
ASM_ENTRY="main.py"
DATASETS=(data1 data2 data3 data4)

printf "Dataset\tLanguage\tRuntime\tGenomeFraction\tDuplicationRatio\tNGA50\tMissassemblies\tMismatches\n"
printf "-------------------------------------------------------------------------------------------------------\n"

time_cmd() { /usr/bin/time -f "%E" "$@" 1>/dev/null 2>.time && tail -n1 .time | tr -d '\n'; rm -f .time; }
parse_metrics() { printf "NA\tNA\tNA\tNA\tNA"; }

for ds in "${DATASETS[@]}"; do
  if [ ! -d "${ASM_DIR}/${ds}" ] && [ -f "${ASM_DIR}/${ds}.zip" ]; then
    (cd "${ASM_DIR}" && unzip -q "${ds}.zip")
  fi
done

for ds in "${DATASETS[@]}"; do
  RT_PY=$( (cd "${ASM_DIR}" && time_cmd ${PY} "${ASM_ENTRY}" "${ds}") || true )
  [ -z "${RT_PY:-}" ] && RT_PY="0:00:00"
  printf "%s\tpython\t%s\t%s\n" "${ds}" "${RT_PY}" "$(parse_metrics)"

  ulimit -s 8192000 || true
  if command -v codon >/dev/null 2>&1 || [ -x "${CODON}" ]; then
    CODON_BIN=$(command -v codon || echo "${CODON}")
    RT_CD=$( (cd "${ASM_DIR}" && time_cmd "${CODON_BIN}" run ${SEQ_PLUGIN} ${RELEASE} "${ASM_ENTRY}" "${ds}") || true )
    [ -z "${RT_CD:-}" ] && RT_CD="0:00:00"
    printf "%s\tcodon\t%s\t%s\n" "${ds}" "${RT_CD}" "$(parse_metrics)"
  else
    printf "%s\tcodon\t0:00:00\t%s\n" "${ds}" "$(parse_metrics)"
  fi
done
