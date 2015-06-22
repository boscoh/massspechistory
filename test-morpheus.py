import logging
from massspechistory import morpheus

logging.basicConfig(level=logging.DEBUG)

morpheus.run({
    '-d': "test-morpheus/ecoli2.mzML",
    '-db': 'db/E_coli_uniprot_iRT.fasta',
    '-ad': 'true',
    '-mmu': 'true',
    '-o': 'test-morpheus/test',
})
