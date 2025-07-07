# BSCAMPP v1.0.8
1. Fixed the path issue with recompiling binaries on macOS.
2. Added `cmake`, `EPA-ng`, and `pplacer` to external requirements when on
macOS systems.
3. Fixed a bug that failed to recover input sequence names when using pplacer
for placements. Each entry in pplacer may contain multiple query sequences,
and previously only the sequence at index 0 is successfully recovered.
3. Now BSCAMPP warns the user if the executables of `EPA-ng` and `pplacer`
are not executable (e.g., `exec format error`), exits the software in
such cases, and notifies the user on how to provide a working executable.

# BSCAMPP v1.0.7
1. Fixed subtree formatting with pplacer output and the final output formatting
regarding the difference between "n" and "nm" in jplace.
2. Added UserWarning suppression to suppress warnings during subtree unifurcation
suppression.

# BSCAMPP v1.0.6
1. Implemented functionality to recompile C++ codes locally so that dynamic
libraries are referred to correctly. 

# BSCAMPP v1.0.5
1. Changed from using original sequence names from query input to using
pseudo names (indexed from 1 to N, N the number of queries). This is to prevent
issues with nontypical characters in sequence names, such as ":".
2. Added better parallelization for the placement phase. Now jobs can run
in parallel to fully take advantage of more cores. By default, each job can
use 2 cores (controlled by `--cpus-per-job`). The number of workers
that can run in parallel is determined by `max(1, num_cpus // cpus_per_job)`.

# BSCAMPP v1.0.4
1. Supported .gz format input files for alignment.
2. Included all other changes made in v1.0.4a.

# BSCAMPP v1.0.4a
1. Removed redundant parameter `--model` since the user-provided info file
is sufficient to run `pplacer` or `epa-ng`. Also this means that it is
not needed to distinguish between nucleotide/protein sequences as long as
the user provides the files correctly.
2. Added a new test for using BSCAMPP on protein sequences, see
`examples/run_bscampp_prot.sh`.

# BSCAMPP v1.0.3
1. Version bump to avoid the entire "pre-release" tag with beta versions.

# BSCAMPP v1.0.2b
1. Removed redundant dependency in `bscampp/jobs.py`.
2. Added logging to each placement subtask with the base method.
3. Changed the temporary file writing directory from a single directory, to
their corresponding directories.

# BSCAMPP v1.0.2
1. Added SCAMPP funtionality and its binary executables.

# BSCAMPP v1.0.1
1. Bumped version to full release.
2. Completed examples for display in `bscampp --help`.

# BSCAMPP v1.0.1b
1. Removed redundant codes and fixed missing variables.
2. Added badges for PyPI installation and current Python Build, etc.
3. Added a single Pytest for dry-running BSCAMPP under `tests/`.
4. For 3, added the `dry_run` parameter to function `pipeline.py:bscampp_pipeline`.

# BSCAMPP v1.0.1a
1. Completed features with both `epa-ng` and `pplacer` support.
2. Refactorized all codes and worked out the PyPI installation for release.
