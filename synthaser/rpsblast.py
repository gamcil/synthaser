import logging
import shutil
import subprocess

from pathlib import Path


LOG = logging.getLogger(__name__)


def get_program_path(program):
    """Get full path to a program on system PATH."""
    path = shutil.which(program)
    if not program:
        raise OSError(f"{program} not found on system $PATH")
    return Path(path).resolve()


def rpsblast(query, database, cpu=2):
    """Run rpsblast on a query file against a database."""
    path = get_program_path("rpsblast")
    params = {
        "-db": database,
        "-comp_based_stats": "1",
        "-seg": "no",
        "-evalue": "3",
        "-num_threads": str(cpu),
        "-outfmt": "11",
    }
    if isinstance(query, str) and Path(query).exists():
        params["-query"] = query
    params = [value for pair in params.items() for value in pair]
    process = subprocess.run(
        [path, *params],
        input=query if "-query" not in params else None,
        stdout=subprocess.PIPE,
    )
    return process


def rpsbproc(results):
    """Convert raw rpsblast results into CD-Search results using rpsbproc.

    Note that since rpsbproc is reliant upon data files that generally are installed in
    the same directory as the executable (and synthaser makes no provisions for them
    being stored elsewhere), we must make sure we have the full path to the original
    executable. If it is called via e.g. symlink, rpsbproc will not find the data files
    it requires and throw an error.

    The CompletedProcess returned by this function contains a standard CD-Search results
    file, able to be parsed directly by the results module.
    """
    path = get_program_path("rpsbproc")
    command = [str(path), "-m", "full", "--quiet", "-t", "doms", "-f"]
    if isinstance(results, str) and Path(results).exists():
        command.extend(["-i", results])
    process = subprocess.run(
        command,
        input=results if "-i" not in command else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    return process


def search(query, database, cpu=2):
    """Convenience function for running rpsblast and rpsbproc."""
    return rpsbproc(rpsblast(query, database, cpu=cpu).stdout)
