# Data

ModelCity does not require bundled data. Pass any suitable pandas DataFrame to
`modelcity.as_settlements()`.

The [`examples/`](examples/) directory contains the two project datasets used
by the convenience readers, integration tests, and demonstration runner:

- `examples/yautepec.csv`: interval records with occupied area and BP dates;
- `examples/viabundus.csv`: snapshot records with population and CE dates.

Treat these files as worked examples rather than a canonical input schema. See
[Adapting a new dataset](../docs/adapting-datasets.md) for role mapping,
scaling, date handling, and validation.

Before redistributing or publishing the example data, add source citations,
licenses, and derivation details here.
