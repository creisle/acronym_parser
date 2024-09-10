# Acronym Parser

A simple regex-based acronym parser to pull of acronym definitions from scientific articles which follow the pattern of

```text
some words definining an acronym (SWDA)
```

## Getting Started

Currently this pacakge must be installed from this repo but in the future I may publish it to pip. I use poetry to install and manage dependencies. Note this requires python3.11 or higher

```bash
git clone ....acronym_parser
cd arconym_parser
poetry install
```

After install you can use the package to parse acronym definitions from bioc documents. In the following example we are downloading a bioc document from pubmed and then applying the acronym parser

```python
import requests

from acronym_parser import mark_acronyms

url = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC8345926/unicode'
resp = requests.get(url)
resp.raise_for_status()

acronyms = mark_acronyms(resp.text)
```
