from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sniffles.config import SnifflesConfig
    from sniffles.sv import SVCall


def get_config(*args, **kwargs) -> 'SnifflesConfig':
    if '--input' not in args:
        args = list(args) + ['--input', 'input.bam']

    if '--vcf' not in args:
        args = list(args) + ['--vcf', 'out.vcf']

    from sniffles.config import SnifflesConfig
    return SnifflesConfig(*args, **kwargs)


def get_svcall(**kwargs) -> 'SVCall':
    args = {
        'contig': 'chr1',
        'pos': 100_000,
        'id': 'Sniffles.Test.1',
        'ref': 'N',
        'alt': 'N',
        'qual': 60,
        'filter': 'PASS',
        'info': {},
        'svtype': 'INS',
        'svlen': 1000,
        'end': 100_001,
        'genotypes': {},
        'precise': True,
        'support': 5,
        'rnames': [],
        'qc': True,
        'nm': -1,
        'postprocess': None,
    }
    args.update(kwargs)

    from sniffles.sv import SVCall
    return SVCall(**args)
