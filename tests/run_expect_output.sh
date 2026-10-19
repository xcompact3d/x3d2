#!/bin/sh
# Usage: run_expect_output.sh <string> <command...>
# Run the command with stderr merged into stdout and print its output.
# Exit with the command status if it is not zero.
# Otherwise exit 1 if the output does not contain the fixed string.
if [ "$#" -lt 2 ]; then
  echo "usage: $0 <string> <command...>" >&2
  exit 2
fi
expected=$1
shift
out=$("$@" 2>&1)
status=$?
printf '%s\n' "$out"
if [ "$status" -ne 0 ]; then
  exit "$status"
fi
printf '%s\n' "$out" | grep -qF -- "$expected"
