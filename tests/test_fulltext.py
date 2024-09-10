import requests
import logging
import pytest
import bioc

from acronym_parser import grab_acronyms

def fetch_biocxml(pmc_id: str, format='unicode') -> str:
    """
    ex. https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC1790863/unicode
    """
    url = (
        f'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmc_id}/{format}'
    )
    logging.info(f'GET {url}')
    resp = requests.get(url)
    resp.raise_for_status()
    return resp.text



@pytest.mark.parametrize('pmc_id, expected_acronyms', [
    pytest.param('PMC8345926', [('ER+', 'Estrogen receptor-positive')], marks=pytest.mark.skip(reason='TODO. support punctuation as part of acronym')),
    ('PMC8345926', [('BC', 'Breast Cancer'), ('ET', 'Endocrine therapy'), ('SERD', 'selective estrogen receptor degrader')])
])
def test_fulltext(pmc_id: str, expected_acronyms: list[tuple[str, str]]):
    content = fetch_biocxml(pmc_id)
    collection = bioc.loads(content)
    doc = collection.documents[0]
    acronyms = grab_acronyms(doc)

    for expected_acronym in expected_acronyms:
        assert expected_acronym in acronyms
