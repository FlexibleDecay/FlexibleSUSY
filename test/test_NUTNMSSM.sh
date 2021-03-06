#!/bin/sh

# directory of this script
BASEDIR=$(dirname $0)

exit_code=0

echo "Running benchmark points BP1-3 from arxiv:1308.1333"

for point in BP1 BP2 BP3
do
    echo "Running $point ..."

    input_file="${BASEDIR}/../model_files/NUTNMSSM/LesHouches.in.NUTNMSSM_1308.1333_${point}"
    output_file="${BASEDIR}/LesHouches.out.NUTNMSSM_${point}"

    ${BASEDIR}/../models/NUTNMSSM/run_NUTNMSSM.x \
        --slha-input-file=${input_file} \
        --slha-output-file=${output_file}

    error="$?"

    expected=0
    # BP2 should fail due to tree-level hh tachyon
    if test "$point" = "BP1" -o "$point" = "BP2"; then
        expected=1
    fi

    if test "x$error" = "x$expected"; then
        echo "=========="
        echo "${point}: OK"
        echo "exit code = $error"
        echo "expected  = $expected"
        echo "=========="
        grep Problems ${output_file}
        grep hh ${output_file}
    else
        echo "=========="
        echo "${point}: FAIL"
        echo "exit code = $error"
        echo "expected  = $expected"
        echo "=========="
        grep Problems ${output_file}
        grep hh ${output_file}
        exit_code=1
    fi

    echo ""

    rm -f ${output_file}
done

exit ${exit_code}
