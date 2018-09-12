#!/bin/bash

# This script installs a clang-format git pre-commit hook for MOOSE, which
# formats all staged files with clang-format prior to a commit, either for a
# stand-alone MOOSE git repository, or for MOOSE as a submodule of another git
# repository. This script can be run from any directory. If the hook file
# already exists, an error message will be printed.

script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# determine the appropriate hook file
function moose_is_git_submodule()
{
    # get the top level of the MOOSE repository or submodule
    local moose_git_dir
    moose_git_dir="$(cd $script_dir && git rev-parse --show-toplevel)"

    # Try running git command one level above the top level MOOSE directory:
    # * If MOOSE is a submodule, then ".." corresponds to the parent repository,
    #   so its top level will be printed.
    # * If MOOSE is NOT a submodule, then ".." leads to a non-git directory, so
    #   git will print an error message that gets redirected to /dev/null, so
    #   nothing will be printed.
    cd $moose_git_dir/.. && git rev-parse --show-toplevel 2> /dev/null
}
if [[ $(moose_is_git_submodule) ]]; then
    hookfile="$(cd $script_dir && git rev-parse --show-toplevel)/../.git/modules/moose/hooks/pre-commit"
else
    hookfile="$(cd $script_dir && git rev-parse --show-toplevel)/.git/hooks/pre-commit"
fi

# ensure that hook file does not exist yet
if [[ -f $hookfile ]]; then
    echo "ERROR: The pre-commit hook '$hookfile' already exists." 1>&2
    echo "       If you would like to replace it, then remove it and run this script again." 1>&2
    exit 1
fi

# create the hook file and make it executable
echo '#!/bin/bash

# Pre-commit hook to run clang-format on all staged c++ files.

patch=$(git clang-format --diff -- $(git diff --staged --name-only -- framework/src framework/include modules/*/src modules/*/include test unit examples tutorials stork))
if [[ "$patch" =~ "no modified files to format" || "$patch" =~ "clang-format did not modify any files" ]]; then
    echo "" > /dev/null
else
    echo ""
    echo "Your code is not properly formatted." >&2
    echo "Run 'git clang-format' to resolve the following issues:" >&2
    echo ""
    echo "$patch"
    exit 1
fi

./scripts/git-hit --staged --check
if [[ $? -ne 0 ]]; then
    echo "One or more HIT files are not properly formatted." >&2
    echo "Run './scripts/git-hit' to resolve the issues." >&2
    exit 1
fi
' > $hookfile

chmod a+x $hookfile
