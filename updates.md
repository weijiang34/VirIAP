# Update logs  

## **v0.2.0** - OVU (Operational Viral Units) updates!  

**Now featuring:** Construction, classification, and abundance estimation of OVUs.
<details>
  <summary>More details</summary>  

**What’s new in v2?**  
This release focuses on **viral contig processing**, with enhanced:  
- ✔ **Clustering** → Group contigs into OVUs  
- ✔ **Mapping** → Assign reads for abundance estimation  
- ✔ **Classification** → Improved taxonomic annotation  

**Key benefits:**
- More accurate viral profiling
- Streamlined workflow for OVU-based analysis
- Better visualization of results
</details>

## **v0.1.0** - Virus Identification and Analysis Pipeline.

**A modular pipeline for virus-specific analysis**, including:  
- Viral signal extraction
- Decontamination (rRNAs from bacteria, archaea, eukaryotes, and mitochondria)
- Deduplication
- Quality assessment & filtering  

<details>
  <summary>More details</summary>  

  ### Modules
  `search`  
  Generates viral detection jobs using four high-performance tools:  
  - CAT_pack
  - VirSorter2
  - GeNomad
  - ViraLM  

  `extract`  
  Selects putative viral contigs based on user-defined criteria.

  `decontam`  
  Removes contaminating rRNA sequences derived from:
  - Bacteria (bac)  
  - Archaea (arc)  
  - Eukaryotes (euk)  
  - Mitochondria (mito)  

  `merge`  
  Combines viral contigs from all samples for downstream processing.

  `dedup`  
  Eliminates duplicate contigs to optimize computational efficiency.

  `check_quality`  
  Filters low-quality viral contigs using quality assessment metrics.

  A summary:  
  | Module         | Functionality                                                                 |  
  |----------------|------------------------------------------------------------------------------|  
  | `search`       | Run viral detection (CAT_pack, VirSorter2, GeNomad, ViraLM)                 |  
  | `extract`      | Select putative viral contigs (user-defined thresholds)                      |  
  | `decontam`     | Filter rRNA contaminants (bac, arc, euk, mito)                              |  
  | `merge`        | Aggregate contigs across samples                                            |  
  | `dedup`        | Remove duplicate contigs                                                    |  
  | `check_quality`| Filter low-quality viral contigs                                            |  
  
</details>  
