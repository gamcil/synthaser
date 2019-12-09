"""
Wrap RPS-BLAST to enable local searches.

Cdd_LE.tar.gz = pre-formatted database containing all domains found in the CDD.

This module requires the following 2 programs:

1) rpsblast
    RPSBLAST is a distributed in the NCBI BLAST+ toolkit. This can be acquired
    either directly from NCBI's FTP:
    ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

    Or from your distributions repositories, e.g. for Ubuntu:
    $ sudo apt install ncbi-blast+

2) rpsbproc
    This is a post-processing tool offered by the NCBI that can convert rpsblast
    results into output identical to that of an online batch CD-Search run.

    This can be installed as follows:
    a) Acquire the relevant archive for your system from the CDD FTP:
        ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/
    b) Extract the contents
    c) Acquire the data files required by rpsbproc either by running
        utils/getcdddata.sh, or directly from the FTP as detailed by the README
        (see: domain-annotation files). The program does NOT require you to
        download all of the domain databases. So, if doing the former, you can
        cancel the run after the necessary files are in data/, then delete db/
        and the database .tar.gz files.
    d) Symlink the rpsbproc binary file to somewhere on your system $PATH.
        This is a requirement of synthaser, as it will throw an error if it
        cannot find rpsbproc directly on the $PATH (i.e. accessable in terminal
        just by typing 'rpsbproc'). For example:
        $ ln -ls /path/to/RpsbProc-x64-linux/rpsbproc /usr/local/bin/rpsbproc
"""


import logging
import shutil
import subprocess
import tarfile

from ftplib import FTP
from pathlib import Path
from datetime import datetime

LOG = logging.getLogger(__name__)

for dependency in ("rpsblast", "rpsbproc"):
    if not shutil.which(dependency):
        raise RuntimeError(f"{dependency} not found on system $PATH")


def download_database(directory, flavour="Cdd"):
    """Download a pre-formatted CDD file from NCBI.

    Database types:
        Cdd:
            All possible domain families
        Cdd_NCBI:
            NCBI-curated families
        Cog:
            Automatically aligned sequences/fragments classified in the COG resource
            (mostly prokaryotes)
        Kog:
            Like Cog, except eukaryote focus
        Pfam:
            Families from the Pfam-A seed alignment database
        Prk:
            Automatically aligned sequences/fragments from the NCBI's Protein Clusters database
        Smart:
            Families from the Smart domain alignment database
        Tigr:
            Families from the TIGRFAM database
    """
    if not Path(directory).is_dir():
        raise TypeError("Invalid directory")

    valid = {"Cdd", "Cdd_NCBI", "Cog", "Kog", "Pfam", "Prk", "Smart", "Tigr"}

    if flavour not in valid:
        raise ValueError("Invalid database flavour")

    db_name = f"{flavour}_LE.tar.gz"
    db_path = Path(directory) / db_name

    LOG.info("Connecting to NCBI FTP server")

    with FTP("ftp.ncbi.nih.gov") as ftp, db_path.open("wb") as db:
        ftp.login()
        ftp.cwd("pub/mmdb/cdd/little_endian/")

        size, date = get_file_facts(db_name, ftp)

        def callback(chunk):
            db.write(chunk)
            progress = db.tell() / size
            print(f"Progress: {progress:.2%}", end="\r", flush=True)

        LOG.info("Attempting to download %s", db_name)
        LOG.info("Size: %.2f gb", size / 1024 ** 3)
        LOG.info("Date: %s", date)

        ftp.retrbinary(f"RETR {db_name}", callback)

    if size != db_path.stat().st_size:
        LOG.warning("Size mismatch between FTP and downloaded copy")

    return db_path


def get_file_facts(database, ftp):
    """Get size and last modified date of a file in FTP directory."""
    for entry, facts in ftp.mlsd(facts=["modify", "size"]):
        if entry == database:
            size = int(facts["size"])
            date = datetime.strptime(facts["modify"], "%Y%m%d%H%M%S").strftime(
                "%d %b %Y %H:%M:%S"
            )
            return size, date
    raise ValueError("Specified database not found in FTP directory")


def untar(filename, dest=None):
    """Untar a file to a specified destination.

    If dest is not specified, will default to using the name of the file without file
    extensions.

    Returns
    -------
    dest: pathlib.Path
        Folder where the contents of the tar file were extracted

    Raises
    ------
    FileNotFoundError
        If file specified by filename does not exist
    """
    filename = Path(filename)

    if not filename.exists():
        raise FileNotFoundError("Given file does not exist")

    if not dest:
        dest = filename.with_suffix("").with_suffix("")
    else:
        dest = Path(dest)

    if not dest.is_dir():
        dest.mkdir()

    with tarfile.open(filename, "r") as tar:
        tar.extractall(dest)

    return dest


def getdb(database, folder):
    """Convenience function to download a collection of databases to a given folder."""
    db_path = download_database(folder, flavour=database)

    LOG.info("Extracting contents to: %s", db_path.with_suffix("").with_suffix(""))
    untar(db_path)


def rpsblast(query, database, cpu=2):
    """Run rpsblast on a query file against a database."""
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
        ["rpsblast", *params],
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
    file, able to be parsed directly by the `results` module.
    """
    path = Path(shutil.which("rpsbproc")).resolve()

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
