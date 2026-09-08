#!/usr/bin/bash

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

clang-format-21 -i $(find $script_dir/srcBeren/ \( -name "*.cpp" -or -name "*.h" \) -and -not -path "*3rd_party*")
