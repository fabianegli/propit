# Proteomics Pipeline Toolkit

![PyPI version](https://img.shields.io/pypi/v/propit.svg)

Glue and Grease for Proteomics Pipelines.

* PyPI package: https://pypi.org/project/propit/
* Free software: GPLv3
* Documentation: https://propit.readthedocs.io.


```
propit --help
usage: propit [-h] {percolator2flashlfq,sage2flashlfq,flashlfq2pb,comet2perc} ...

Convert Proteomics files.

positional arguments:
  {percolator2flashlfq,sage2flashlfq,flashlfq2pb,comet2perc}
    percolator2flashlfq
                        Make generic input files for FlashLFQ.
    sage2flashlfq       Make generic input files for FlashLFQ.
    flashlfq2pb         Make generic input files for ProteoBench from FlashLFQ peptides.
    comet2perc          This tool removes the deltaCn column from pin files.

options:
  -h, --help            show this help message and exit
```


## Credits

This package was created with [Cookiecutter](https://github.com/audreyfeldroy/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
