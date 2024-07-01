from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sniffles.config import SnifflesConfig


def get_config(*args, **kwargs) -> 'SnifflesConfig':
    if '--input' not in args:
        args = list(args) + ['--input', 'input.bam']

    if '--vcf' not in args:
        args = list(args) + ['--vcf', 'out.vcf']

    from sniffles.config import SnifflesConfig
    return SnifflesConfig(*args, **kwargs)