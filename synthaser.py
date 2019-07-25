#!/usr/bin/env python3
"""
Parse CD-search results to extract domain architectures of synthases.

Can also extract domain sequences for each synthase (--extract)
and generate a SVG figure of all synthases and their domain architecture
(--visual).

Author: Cameron Gilchrist
Date: 2018-06-12
"""

import click

from collections import defaultdict, namedtuple, Counter


def group_overlapping_hits(hits):
    """ Iterator that groups Hit namedtuples based on
        overlapping locations.
    """
    def overlapping(one, two):
        """ Return True if Hits overlap, checking both directions.
        """
        one_len = one.end - one.start
        two_len = two.end - two.start
        smallest = one if one_len <= two_len else two

        one_set = set(range(one.start, one.end))
        two_set = set(range(two.start, two.end))
        intersect = one_set & two_set

        overlap = len(intersect) / (smallest.end - smallest.start)
        return True if overlap > 0.9 else False

    hits.sort(key=lambda x: x.start)

    i, total = 0, len(hits)
    while i < total:
        current = hits[i]  # grab current hit
        group = [current]  # start group

        if i == total - 1:  # if current hit is the last, yield
            yield group
            break

        for j in range(i + 1, total):  # iterate rest
            future = hits[j]           # grab next hit

            if overlapping(current, future):
                group.append(future)   # add if contained
            else:
                yield group            # else yield to iterator
                break

            if j == total - 1:  # if reached the end, yield
                yield group

        i += len(group)  # move index ahead of last group


def parse_results_file(results):
    """ Extract all hits for every query.

        Args:
            cd_search (str): CD-search results
        Returns:
            domains (dict): Query:[hits]
    """
    domains = {
            'KS': ['PKS_KS', 'PKS'],
            'AT': ['PKS_AT', 'Acyl_transf_1'],
            'ER': ['PKS_ER', 'enoyl_red'],
            'KR': ['KR', 'PKS_KR'],
            'TE': ['Thioesterase', 'Aes'],
            'TR': ['Thioester-redct', 'SDR_e1'],
            'MT': ['Methyltransf_11',
                   'Methyltransf_12',
                   'Methyltransf_23',
                   'Methyltransf_25',
                   'Methyltransf_31',
                   'AdoMet_MTases'],
            'DH': ['PKS_DH', 'PS-DH'],
            'PT': ['PT_fungal_PKS'],
            'ACP': ['PKS_PP', 'PP-binding', 'AcpP'],
            'SAT': ['SAT'],
            'C': ['Condensation'],
            'A': ['A_NRPS', 'AMP-binding']
            }

    Hit = namedtuple('Hit', ['type', 'domain', 'start', 'end'])
    queries = defaultdict(list)
    for row in results:
        row = row.decode()
        if not row.startswith('Q#'):
            continue
        else:
            fields = row.split('\t')

        # Query protein ID stored as Q#1 - >proteinID
        query, name = fields[0].split('>')[1], fields[8]
        for domain_type, domain_names in domains.items():
            if name not in domain_names:
                continue
            else:
                start, end = int(fields[3]), int(fields[4])
                hit = Hit(domain_type, name, start, end)
                queries[query].append(hit)
                break
    return queries


def find_architecture(hits):
    """ Determine domain architecture given a list of domains.

        Iterates overlapping domain groups, then saves the
        largest hit of each.
    """
    # Save final domain list for this query
    domains = []
    for group in group_overlapping_hits(hits):

        # Save the largest in each group
        largest = max(group, key=lambda x: x.end - x.start)
        domains.append(largest)

    return domains


def parse_fasta(fasta):
    """ Parse an open FASTA file for sequences.
    """
    sequences = {}
    for line in fasta:
        line = line.strip()
        if line.startswith('>'):
            header = line[1:]
            sequences[header] = ''
        else:
            sequences[header] += line
    return sequences


def generate_SVG(synthases, sequences, spacing=40, width=600):
    """ Build an SVG figure given parsed synthases and sequences.
    """
    def create_linear_gradient(domains, synthase, sequence_length):
        """ Create a coloured gradient based on domain architecture.
        """
        colours = {'ACP': '#084BC6', 'KS': '#08B208',
                   'SAT': '#808080', 'KR': '#089E4B',
                   'MT': '#00ff00', 'ER': '#089E85',
                   'AT': '#DC0404', 'DH': '#B45F04',
                   'PT': '#999900', 'TE': '#750072',
                   'TR': '#9933ff', 'T': '#084BC6',
                   'R': '#9933ff', 'C': '#393989',
                   'A': '#56157F'}

        # Add first stop in the gradient; gap from start of protein
        # to the first domain
        stops = ['<stop offset="0%" stop-color="white" />']

        for domain in domains:
            # Get start and end as percentages of whole
            start_pct = int(domain.start / sequence_length * 100)
            end_pct = int(domain.end / sequence_length * 100)

            # Create gradient stops. Have to do two stops per coordinate
            # (white/colour at start, then colour/white at end), so these
            # will be drawn as hard edges rather than actual gradients
            colour = colours[domain.type]
            stops.append(
                    f'<stop offset="{start_pct}%" stop-color="white" />\n'
                    f'<stop offset="{start_pct}%" stop-color="{colour}" />\n'
                    f'<stop offset="{end_pct}%" stop-color="{colour}" />\n'
                    f'<stop offset="{end_pct}%" stop-color="white" />')

        return ('<linearGradient id="{}_doms" x1="0%" y1="0%" x2="100%"'
                ' y2="0%">\n{}'
                '\n</linearGradient>'''.format(synthase, '\n'.join(stops)))

    def synthase_svg(synthase, domains, sequence, scale_factor):
        """ Build SVG representation of one synthase.

            A----------B
            |           \
            |            C
            |           /
            E----------D
        """
        # Scale sequence length
        sequence_length = len(sequence)
        scaled_length = scale_factor * sequence_length

        # Get coordinates for this synthase
        ax, ay = 0, 5
        bx, by = scaled_length - 10, 5
        cx, cy = scaled_length, 12
        dx, dy = scaled_length - 10, 19
        ex, ey = 0, 19

        # Concatenate for SVG
        points = f'{ax},{ay},{bx},{by},{cx},{cy},{dx},{dy},{ex},{ey}'

        # Create linearGradient based on domain architecture
        gradient = create_linear_gradient(domains, synthase, sequence_length)

        # Form information string to print above the polygon
        architecture = '-'.join(d.type for d in domains)
        information = f'{synthase}, {sequence_length}aa, {architecture}'

        # Create the polygon
        polygon = (f'<text y="0" font-size="12">{information}</text>'
                   f'<polygon id="{synthase}" points="{points}"'
                   f' fill="url(#{synthase}_doms)" stroke="black"'
                   ' stroke-width="1.5"/>')

        return gradient, polygon

    image = ''
    gradients, polygons = [], []

    # Find largest synthase to determine scaling factor
    largest_sequence = max(sequences.values(), key=len)
    scale_factor = (width - 2) / len(largest_sequence)

    for synthase_type in synthases:
        type_polygons = []

        # Descending sort by sequence length
        synthases[synthase_type].sort(key=lambda x: len(sequences[x[0]]),
                                      reverse=True)

        for synthase, domains in synthases[synthase_type]:
            sequence = sequences[synthase]

            # Form the polygon and its fill gradient
            gradient, polygon = synthase_svg(synthase, domains,
                                             sequence,
                                             scale_factor)
            gradients.append(gradient)
            type_polygons.append(polygon)

        polygons.append([synthase_type] + type_polygons)

    # Gradient definitions
    gradients = '\n'.join(gradients)
    image += f'\n<defs>\n{gradients}</defs>'

    # Filter out empty blocks
    polygons = [p for p in polygons if len(p) > 1]
    block_offset = 12
    for type_index, polygons in enumerate(polygons):

        image += (f'\n<g>\n<text x="-1" y="{block_offset}" font-size="15"'
                  f' font-weight="bold">{polygons[0]}</text>')

        for index, polygon in enumerate(polygons[1:]):
            poly_offset = index * spacing
            offset = block_offset + poly_offset + 20
            image += (f'\n<g transform="translate(1,{offset})">'
                      f'\n{polygon}\n</g>')

        # Since we use this for sizing the SVG canvas, check
        # if it's the final block
        if type_index != len(polygons) - 1:
            block_offset += poly_offset + 2 * spacing

        image += '\n</g>'

    # End tag for svg and return
    block_offset -= 40
    header = f'<svg width="{width}" height="{block_offset}">'
    footer = '\n</svg>'
    return header + image + footer


def wrap_fasta(sequence, limit=80):
    """ Wrap FASTA record to 80 characters per line.
    """
    return '\n'.join(sequence[i: i + limit] for
                     i in range(0, len(sequence), limit))


def classify_synthases(synthase_type, synthases):
    """ Sort synthases into specific subtypes based on domain composition.
    """
    if synthase_type == 'pks':
        classified = {'HR-PKS': [], 'NR-PKS': [], 'PR-PKS': [], 'Other': []}
        required = {
            'HR-PKS': {'ER', 'KR', 'DH'},
            'PR-PKS': {'KR', 'DH'},
            'NR-PKS': {'SAT', 'PT'},
            'Other': set()
        }
        for query, domains in synthases.items():
            combo = (query, domains)
            types = set(d.type for d in domains)
            for pks in ('HR-PKS', 'PR-PKS', 'NR-PKS', 'Other'):
                if required[pks].issubset(types):
                    classified[pks].append(combo)
                    continue

    elif synthase_type == 'nrps':
        classified = {'NRPS': [], 'NRPS-like': [], 'Other': []}
        required = {
            'NRPS': {'A', 'T', 'C'},
            'NRPS-like': {'A'},
            'Other': set()
        }
        for query, domains in synthases.items():
            # Replace ACP with T, TR with R as is convention
            for index, domain in enumerate(domains):
                if domain.type == 'ACP':
                    domains[index] = domain._replace(type='T')
                elif domain.type == 'TR':
                    domains[index] = domain._replace(type='R')

            combo = (query, domains)
            types = set(d.type for d in domains)
            for nrps in ('NRPS', 'NRPS-like'):
                if required[nrps].issubset(types):
                    classified[nrps].append(combo)
                    continue

    elif synthase_type == 'hybrid':
        classified = {'PKS-NRPS': []}
        for query, domains in synthases.items():
            condensation = False
            for index, domain in enumerate(domains):
                # Condensation domain marks the start of the NRPS module
                if not condensation:
                    if domain.type == 'C':
                        condensation = True
                    else:
                        continue
                # Replace ACP with T, TR with R, as above
                if domain.type == 'ACP':
                    domains[index] = domain._replace(type='T')
                elif domain.type == 'TR':
                    domains[index] = domain._replace(type='R')
            combo = (query, domains)
            classified['PKS-NRPS'].append(combo)

    return classified


@click.command()
@click.argument('synthase_type',
                type=click.Choice(['pks', 'nrps', 'hybrid']))
@click.argument('results', type=click.File('rb'))
@click.option('-f', '--fasta',
              help='Proteins submitted to CD-search (FASTA)',
              type=click.File('rb'))
@click.option('-v', '--visualise',
              help='Generate .svg visualisation',
              default=False)
@click.option('-e', '--extract',
              help='Extract domain regions from each synthase',
              default=False)
@click.option('-o', '--output',
              help='Base name for output files')
def synthaser(synthase_type, results, fasta, visualise, extract, output):
    """ Parse NCBI CD-Search results file to find domain architectures
        of secondary metabolite synthases.
    """
    # Parse CD search results file for domain hits
    queries = parse_results_file(results)

    # Filter out overlaps, find domain architecture
    for query, domains in queries.items():
        queries[query] = find_architecture(domains)

    # Classify synthases based on domain composition
    classified = classify_synthases(synthase_type, queries)

    # Print out domain architectures by synthase type
    for i, synthase_type in enumerate(classified):
        print(synthase_type)
        print('-' * len(synthase_type))
        for synthase, domains in classified[synthase_type]:
            print(synthase, ' ', '-'.join(x.type for x in domains))
        if i == 3:
            break
        print('')

    if fasta:
        # Parse FASTA file, build dictionary mapping sequences
        # to corresponding protein IDs
        sequences = parse_fasta(fasta)

    if visualise:
        # Generate text of an SVG figure
        image = generate_SVG(classified, sequences)

        # Write to file
        with open(f'{output}.svg', 'w') as output:
            output.write(image)

    if extract:
        with open(f'{output}_domains.faa', 'w') as out:
            for synthase_type in classified:
                for synthase, domains in classified[synthase_type]:
                    counter = Counter()
                    for domain in domains:
                        # Add to counter, then use the current count
                        # to dynamically form a numbered header
                        counter.update(domain.type)
                        header = '>{}_{}_{}'.format(synthase,
                                                    domain.type,
                                                    counter[domain.type])

                        # Slice domain sequence
                        sequence = sequences[synthase][domain.start-1:
                                                       domain.end]

                        # Format wrapped FASTA and write
                        fasta = '{}\n{}\n'.format(header,
                                                  wrap_fasta(sequence))
                        out.write(fasta)


if __name__ == '__main__':
    synthaser()
