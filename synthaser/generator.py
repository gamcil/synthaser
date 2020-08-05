"""Domain rule file generator.

"""

import argparse
import json
import gzip
import shutil

from ftplib import FTP
from pathlib import Path

from synthaser.rpsblast import download, untar


def get_parser():
    parser = argparse.ArgumentParser("rulegen")
    parser.add_argument(
        "domain_file",
        type=argparse.FileType("r+"),
        help="Domain rule file",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--list",
        help="List all current domain rules",
        action="store_true",
    )
    group.add_argument(
        "--print",
        help="Print given domain entry",
        nargs="+"
    )
    group.add_argument(
        "--update",
        type=argparse.FileType("r"),
        help="Add new or update existing domain rules",
    )
    group.add_argument(
        "--delete",
        help="Delete rules",
        nargs="+",
    )
    group.add_argument(
        "--download",
        help="Download CDD family files",
    )
    group.add_argument(
        "--parse",
        help="Folder containing CDD family files"
    )
    return parser


def list_domains(d):
    offset = max(len(c["type"]) for c in d["classes"])
    for item in d["classes"]:
        if "families" in item:
            fams = len(item["families"])
        else:
            fams = 0
        line = "{}  {} [{} families]".format(
            item["type"].rjust(offset),
            item["name"],
            fams
        )
        print(line)


def print_domains(d, domains):
    for item in d["classes"]:
        if item["type"] in domains:
            x = json.dumps(item, indent=2)
            print(x)


def delete_domains(d, domains):
    d["classes"] = [
        domain
        for domain in d["classes"]
        if domain["type"] not in domains
    ]


def parse_domains(fp):
    entries = []
    for entry in fp.read().split(">"):
        if not entry:
            continue
        header, *families = entry.strip().split("\n")
        acronym, name = header.split(" ", 1)
        d = dict(type=acronym, name=name, families=families)
        entries.append(d)
    return entries


def download_families(entries):
    families = []
    for entry in entries:
        families.extend(entry["families"])

    print("families", families)



def download_cdd_families(folder):
    folder = Path(folder)

    cddid_gz = folder / "cddid_all.tbl.gz"
    cddid = cddid_gz.with_suffix("")

    bitscore = folder / "bitscore_specific.txt"
    families = folder / "family_superfamily_links"

    with FTP("ftp.ncbi.nih.gov") as ftp:
        ftp.login()
        ftp.cwd("pub/mmdb/cdd")

        for path in [cddid_gz, bitscore, families]:
            print(ftp, path.name, path)
            download(ftp, path.name, path)

        with gzip.open(cddid_gz, "rb") as f_in, cddid.open("wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        cddid_gz.unlink()


def parse_cdd_families(folder, indent=2):
    """Generates CDD domain family information database.

    Parses three files downloaded from the NCBI FTP:
    1. cddid_all.tbl: Family information, PSSM length
    2. bitscore_specific.txt: Bitscore thresholds for specific hits
    3. family_superfamily_links: Family-Superfamily relationships

    Then dumps to JSON. File can be used when updating the domain rule
    JSON file with this module.
    """
    folder = Path(folder)

    cddid = folder / "cddid_all.tbl"
    bitscore = folder / "bitscore_specific.txt"
    families = folder / "family_superfamily_links"
    summary = folder / "summary.json"

    d = {}

    with cddid.open() as fp:
        parse_cddid(fp, d)

    with bitscore.open() as fp:
        parse_bitscores(fp, d)

    with families.open() as fp:
        parse_families(fp, d)

    with summary.open("w") as fp:
        json.dump(d, fp, indent=indent)


def parse_cddid(fp, d):
    for line in fp:
        if not line or line.isspace():
            continue
        pssm, accession, name, _, length = line.strip().split("\t")
        d[accession] = dict(
            pssm=pssm,
            accession=accession,
            name=name,
            length=length,
        )


def parse_bitscores(fp, d):
    for line in fp:
        if not line or line.isspace():
            continue
        _, accession, bitscore = line.strip().split("\t")
        d[accession]["bitscore"] = bitscore


def parse_families(fp, d):
    for line in fp:
        if not line or line.isspace():
            continue
        family, _, superfamily, _ = line.strip().split("\t")
        d[family]["superfamily"] = superfamily


def update_domains(d, domains):
    """Adds new, or updates existing domain rules."""
    entries = parse_domains(domains)
    download_families(entries)
    pass


def main(args):
    d = json.load(args.domain_file)

    if args.list:
        list_domains(d)
        return

    if args.print:
        print_domains(d, args.print)
        return

    if args.download:
        download_cdd_families(args.download)
        return

    if args.parse:
        parse_cdd_families(args.parse)
        return

    if args.delete:
        delete_domains(d, args.delete)

    elif args.update:
        update_domains(d, args.update)

    print(f"Writing: {args.domain_file.name}")
    # json.dump(d, args.domain_file)


if __name__ == "__main__":
    p = get_parser()
    args = p.parse_args()
    main(args)
