#!/usr/bin/env bash

MAIN_BRANCH=${1:-"origin/main"}
CLANG_BIN=${2:-"`which clang-format-14`"}
COMMAND_OUTPUT=$(git-clang-format --diff $MAIN_BRANCH --extensions c,cpp,h,hpp --binary $CLANG_BIN)
grep '^no modified files to format$' <<<"$COMMAND_OUTPUT" && exit 0
grep '^clang-format did not modify any files$' <<<"$COMMAND_OUTPUT" && exit 0
echo -e "\n/// code format check failed! \n"
echo "$COMMAND_OUTPUT"
exit 1
