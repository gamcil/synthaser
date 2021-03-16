"""Extract CDS sequences from GenBank files for synthaser analysis."""

import logging
import re

from pathlib import Path

from Bio import SeqIO


LOG = logging.getLogger(__name__)


def get_feature_type(feature):
    """Tests if an antiSMASH feature is a NRPS or PKS.

    Looks for a /NRPS_PKS="type:..." field containing either 'PKS' or 'NRPS'.
    This function returns True if a matching qualifier is found.
    """
    try:
        ftype = re.search(
            r"type: ([A-Za-z-_\s]+?)$",
            feature.annotations["NRPS_PKS"],
            re.MULTILINE
        ).group(1)
    except (KeyError, IndexError, AttributeError):
        return None
    if "PKS" in ftype:
        return "PKS"
    if "NRPS" in ftype:
        return "NRPS"


def get_NRPS_PKS(features):
    pks, nrps = [], []
    for feature in features:
        ftype = get_feature_type(feature)
        if ftype == "PKS":
            pks.append(feature)
        elif ftype == "NRPS":
            nrps.append(ftype)
    return pks, nrps


def write(path, features):
    with open(path, "w") as fp:
        LOG.info("Writing %i feature(s) to file: %s", len(features), fp.name)
        SeqIO.write(features, fp, "fasta")


def convert(path, antismash=False):
    """Extracts proteins from a GenBank file.

    Arguments:
        path (str): GenBank file path
        antismash (bool): Only save PKS/NRPS sequences
    """
    LOG.info("Parsing GenBank file: %s", path)

    path = Path(path)

    # Parse CDS features from the file
    with path.open() as fp:
        features = [record for record in SeqIO.parse(fp, "genbank-cds")]

    # Blank out descriptions for clean FASTA headers
    for feature in features:
        feature.description = ""

    # If antismash=True, look for PKS and NRPS sequences only
    if antismash:
        LOG.info("Finding antiSMASH PKS and NRPS features")
        pks, nrps = get_NRPS_PKS(features)
        for fts, text in [(pks, 'pks'), (nrps, 'nrps')]:
            if not fts:
                LOG.info("No %s features found, skipping", text.upper())
                continue
            write(path.with_name(f"{path.stem}_{text}.fa"), fts)
    else:
        write(path.with_suffix(".fa"), features)
