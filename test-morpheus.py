import logging
from massspechistory import morpheus

logging.basicConfig(level=logging.DEBUG)

morpheus.run({
    '-d': "test-morpheus/ecoli.raw",
    '-db': 'db/E_coli_uniprot_iRT.fasta',
    '-ad': 'true',
    '-mmu': 'true',
    '-o': 'test-morpheus/test',
    '-vm': 'Ox',
    '-fm': 'AlkC',
    '-acs': 'false'
})

