# 043 — Millisecond Encoding Precision in Auditory Cortex

## Metadata
- **Title:** Millisecond Encoding Precision of Auditory Cortex Neurons
- **Authors:** Christoph Kayser, Nikos K Logothetis, Stefano Panzeri
- **Year:** 2010
- **Journal:** Proceedings of the National Academy of Sciences (PNAS)
- **DOI:** 10.1073/pnas.1012656107
- **PubMed:** 20837521
- **URL:** https://www.pnas.org/doi/10.1073/pnas.1012656107

## Abstract
The research examined how neurons in the auditory cortex of alert primates encode sound information with precise spike timing. Researchers found that registering spikes at a precision coarser than a few milliseconds significantly reduced the encoded information. The study demonstrates that precise neural firing patterns carry stimulus information about complex sounds, with the information loss depending on temporal precision. Rapid firing rate changes — rather than complex spike correlations — were identified as the primary mechanism for fine-timed encoding.

## Key Findings
- Auditory cortex neurons encode information with millisecond precision
- Degrading temporal precision to >few milliseconds significantly reduces encoded information
- Rapid firing rate changes (not spike correlations) drive the fine-timed encoding
- Millisecond precision is fundamental throughout auditory processing systems
- Complex natural sounds require high temporal precision for full information extraction

## Relevance to Aureonoise
The auditory cortex operates at millisecond precision. This means Aureonoise's temporal parameters (grain onset, bilateral alternation timing, envelope attack) must be controlled at the millisecond level to fully exploit the auditory pathway. Coarser temporal control would lose information. This validates the research finding in MEMORY.md: bilateral stimulation requires onset <5ms. The DSP engine must maintain sample-accurate timing — any jitter >1ms degrades the neural encoding of the stimulus.
