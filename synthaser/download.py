"""Database download module."""

import json
import gzip
import shutil
import tempfile
import logging
import tarfile

from ftplib import FTP
from pathlib import Path
from datetime import datetime


LOG = logging.getLogger(__name__)


def download_from_ftp(ftp, retr, out_path):
    with out_path.open("wb") as out:
        size, date = get_file_facts(retr, ftp)

        def callback(chunk):
            out.write(chunk)
            progress = out.tell() / size
            print(f"Progress: {progress:.2%}", end="\r", flush=True)

        LOG.info("Attempting to download %s", retr)
        LOG.info("Size: %.2f gb", size / 1024 ** 3)
        LOG.info("Date: %s", date)

        ftp.retrbinary(f"RETR {retr}", callback)

    if size != out_path.stat().st_size:
        LOG.warning("Size mismatch between FTP and downloaded copy")


def download_database(directory, flavour="Cdd"):
    """Download a pre-formatted CDD file from NCBI.

    Cdd: All possible domain families
    Cdd_NCBI: NCBI-curated families
    Cog: Automatically aligned sequences/fragments classified in the COG resource
        (mostly prokaryotes)
    Kog: Like Cog, except eukaryote focus
    Pfam: Families from the Pfam-A seed alignment database
    Prk: Automatically aligned sequences/fragments from the NCBI's Protein Clusters database
    Smart: Families from the Smart domain alignment database
    Tigr: Families from the TIGRFAM database
    """
    if not Path(directory).is_dir():
        raise TypeError("Invalid directory")

    valid = {"Cdd", "Cdd_NCBI", "Cog", "Kog", "Pfam", "Prk", "Smart", "Tigr"}

    if flavour not in valid:
        raise ValueError("Invalid database flavour")

    db_name = f"{flavour}_LE.tar.gz"
    db_path = Path(directory) / db_name

    LOG.info("Connecting to NCBI FTP server")

    with FTP("ftp.ncbi.nih.gov") as ftp:
        ftp.login()
        ftp.cwd("pub/mmdb/cdd/little_endian/")
        download_from_ftp(ftp, db_name, db_path)

    return db_path


def get_file_facts(database, ftp):
    """Get size and last modified date of a file in FTP directory.

    Parameters:
        database (str): Name of CDD database to check
        ftp (str): FTP connection
    Returns:
        size (int): Size of file
        date (str): Last modified date of file
    Raises:
        ValueError: If database is not found in the FTP directory
    """
    for entry, facts in ftp.mlsd(facts=["modify", "size"]):
        if entry == database:
            size = int(facts["size"])
            date = datetime.strptime(facts["modify"], "%Y%m%d%H%M%S").strftime(
                "%d %b %Y %H:%M:%S"
            )
            return size, date
    raise ValueError("Specified file not found in FTP directory")


def untar(filename, dest=None):
    """Untar a file to a specified destination.

    If dest is not specified, will default to using the name of the file without file
    extensions.

    Returns:
        dest (pathlib.Path): Folder where the contents of the tar file were extracted
    Raises:
        FileNotFoundError: If file specified by filename does not exist
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


def download_cdd_family_files(folder):
    """Downloads CDD family information files from NCBI FTP to a folder."""
    folder = Path(folder)

    if not folder.exists() or not folder.is_dir():
        folder.mkdir()

    cddid_gz = folder / "cddid_all.tbl.gz"
    cddid = cddid_gz.with_suffix("")

    bitscore = folder / "bitscore_specific.txt"
    families = folder / "family_superfamily_links"

    LOG.info("Connecting to NCBI FTP")
    with FTP("ftp.ncbi.nih.gov") as ftp:
        ftp.login()
        ftp.cwd("pub/mmdb/cdd")

        for path in [cddid_gz, bitscore, families]:
            download_from_ftp(ftp, path.name, path)

        LOG.info("Extracting cddid.tbl.gz")
        with gzip.open(cddid_gz, "rb") as f_in, cddid.open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        LOG.info("Deleting cddid.tbl.gz")
        cddid_gz.unlink()


def parse_cdd_families(folder):
    """Generates CDD domain family information database.

    Parses three files downloaded from the NCBI FTP:
    1. cddid_all.tbl: Family information, PSSM length
    2. bitscore_specific.txt: Bitscore thresholds for specific hits
    3. family_superfamily_links: Family-Superfamily relationships

    Then dumps to JSON. File can be used when updating the domain rule
    JSON file with this module.
    """
    folder = Path(folder)

    d = {}
    LOG.info("Building CDD family file...")

    with (folder / "cddid_all.tbl").open() as fp:
        parse_cddid(fp, d)
    with (folder / "bitscore_specific.txt").open() as fp:
        parse_bitscores(fp, d)
    with (folder / "family_superfamily_links").open() as fp:
        parse_families(fp, d)

    return d


def parse_cddid(fp, d):
    """Parses cddid_all.tbl."""
    for line in fp:
        if not line or line.isspace():
            continue
        pssm, accession, name, _, length = line.strip().split("\t")
        d[accession] = dict(
            pssm=int(pssm),
            accession=accession,
            name=name,
            length=int(length),
        )


def parse_bitscores(fp, d):
    """Parses bitscore_specific.txt."""
    for line in fp:
        if not line or line.isspace():
            continue
        _, accession, bitscore = line.strip().split("\t")
        d[accession]["bitscore"] = float(bitscore)


def parse_families(fp, d):
    """Parses family_superfamily_links."""
    for line in fp:
        if not line or line.isspace():
            continue
        family, _, superfamily, _ = line.strip().split("\t")
        d[family]["superfamily"] = superfamily


def download_cdd_families(folder=None, indent=2):
    """Download CDD family information."""
    path = Path('cdd_families.txt')
    if folder:
        download_cdd_family_files(folder)
        d = parse_cdd_families(folder)
        path = folder / path
        with path.open("w") as fp:
            json.dump(d, fp, indent=indent)
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            download_cdd_family_files(tmpdir)
            d = parse_cdd_families(tmpdir)
        with open(path, "w") as fp:
            json.dump(d, fp, indent=indent)
    return path


def getdb(database, folder=None):
    """Main function for download module.
    """
    LOG.info("Starting getdb module")
    if database == "cdd_families":
        path = download_cdd_families(folder=folder)
    else:
        path = download_database(folder, flavour=database)
        untar(path)
    LOG.info("Database written to: %s", path)
