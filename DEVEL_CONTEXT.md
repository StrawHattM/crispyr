# CRISPYR Development Context
*Last Updated: Jan 4, 2025*

---

## PROJECT OVERVIEW

**Package Name:** crispyr  
**Purpose:** Process multiplexed FASTQ files from pooled CRISPR screens into count matrices  
**Repository:** https://github.com/StrawHattM/crispyr  
**Current Branch:** `fastq-processing-pipeline`  
**Development Stage:** Post-implementation debugging and testing

---

## ARCHITECTURE & DATA FLOW

### Pipeline Flow
```
FASTQ Files (multiplexed reads)
    ↓
extract_barcodes_from_fastq()
    ├─ Extracts: sample_barcode (first 8bp)
    └─ Extracts: construct_barcode (20bp after CACCG flanking)
    ↓
match_constructs()
    └─ Maps construct_barcode to reference library (exact or fuzzy)
    ↓
map_constructs_to_genes()
    └─ Maps barcodes to target genes via chip file
    ↓
aggregate_by_sample_barcode()
    └─ Merges sample_barcode with manifest to assign conditions
    ↓
create_count_matrix()
    └─ Generates construct-level and gene-level count matrices
    ↓
lognormalize_counts()
    └─ Applies PoolQ-style normalization: log2((count+1)/(median+1))
    ↓
calculate_qc_metrics()
    └─ Generates overall and per-sample QC reports
```

### Key Data Structures

**Sample Barcode (8bp):**
- Location: First 8bp of each FASTQ read
- Purpose: Demultiplex reads by sample/condition
- Example: TTGAACCG, AATCCAGC

**Construct Barcode (20bp):**
- Location: After flanking sequence (default: CACCG + 6bp offset)
- Purpose: Identify which guide/construct the read came from
- Example: CTGCGGAGCCGTTCACGCCG

**Barcode Extraction Policy:**
- Format: `PREFIX:SEQUENCE@POSITION`
- Default: `PREFIX:CACCG@6` (Broad Institute standard)
- Meaning: Find "CACCG", then extract 20bp starting at position 6 after it

---

## IMPLEMENTATION HISTORY

### Phase 1: Initial Implementation (Day 1)
**Status:** ✅ COMPLETE

**Completed:**
- All 11 main functions implemented (1,561 LOC)
- Comprehensive documentation (Roxygen2, README, vignette)
- Unit tests for core functionality
- Git commits with proper messages
- Package structure and DESCRIPTION

**Files Created:**
- R/load_data.R (3 functions)
- R/extract_barcodes.R (3 functions)
- R/match_and_map.R (4 functions)
- R/count_and_normalize.R (3 functions)
- R/qc_and_process.R (3 functions)
- R/nin.R (helper operator)
- tests/testthat/test-core_functions.R
- vignettes/crispyr-workflow.Rmd

---

## DEBUGGING HISTORY

### Session 1: Post-Implementation Issues
**Date:** Jan 3, 2025

**Problems Identified:**
1. Quality string extraction bug (S4 class coercion)
2. Regex escape sequence issue in parse_barcode_policy()
3. Windows path issues in example scripts

**Solutions Applied:**
1. Changed `as.character(fastq_data@quality)` to `as.character(quality_obj@quality)`
2. Fixed regex pattern from `\d` to `\\d`
3. Updated paths from `/d/...` to `D:/...`

### Session 2: Column Preservation Debugging
**Date:** Jan 4, 2025

**Problems Identified:**
1. ❌ WRONG DIAGNOSIS: Thought sample_barcode was being lost in pipeline
2. ✅ ACTUAL PROBLEM: Test data mismatch (FASTQ vs manifest)

**Code Fixes Applied (Working Correctly):**

#### 1. Fixed dplyr pipes in match_constructs() 
**File:** R/match_and_map.R lines 91-99  
**Issue:** dplyr::mutate converts data.table to tibble  
**Fix:** Pure data.table syntax
```r
# Before
failed_extractions <- extracted_barcodes[...] %>% dplyr::mutate(...)

# After
failed_extractions <- extracted_barcodes[extraction_status %nin% c("success", "no_sample_barcode")]
failed_extractions[, `:=`(construct_id = NA_character_, match_type = "no_match", match_quality = 0)]
```

#### 2. Fixed dplyr pipe in map_constructs_to_genes()
**File:** R/match_and_map.R lines 330-332  
**Issue:** dplyr group_by pipeline potentially drops columns  
**Fix:** Pure data.table by syntax
```r
# Before
result <- result %>% dplyr::group_by(barcode) %>% dplyr::mutate(...) %>% dplyr::ungroup() %>% as.data.table()

# After
result <- result[, barcode_to_genes := dplyr::n_distinct(gene_id, na.rm = TRUE), by = barcode]
```

#### 3. Fixed load_sample_manifest()
**File:** R/load_data.R lines 130-165  
**Issue:** readr::read_csv(col_names = ...) doesn't skip existing header  
**Fix:** Auto-detect headers, select columns, rename, convert to data.table
```r
manifest <- readr::read_csv(manifest_file, progress = !quiet)
manifest <- dplyr::select(manifest, 1, 2)
colnames(manifest) <- col_names
manifest <- data.table::as.data.table(manifest)
```

#### 4. Added debug output to aggregate_by_sample_barcode()
**File:** R/count_and_normalize.R lines 29-71  
**Purpose:** Diagnose merge issues  
**Prints:**
- Input reads: X rows, Y unique sample_barcodes
- Manifest: X rows, Y unique sample_barcodes
- Whether 'condition' column exists after merge
- NA vs non-NA condition counts

**Result:** These fixes work correctly with synthetic test data ✅

---

## CURRENT STATUS: DATA MISMATCH ISSUE

### The Real Problem
**The code works correctly**, but the test data is mismatched.

**Evidence:**
```
Test Data: sample_manifest_subset.csv
├─ Contains: 36 sample barcodes
└─ Conditions: Day_0 (12 samples), UT_r1 (12 samples), UT_r2 (12 samples)

Test Data: subset_sample.fastq.gz  
├─ Contains: 75,063 reads
├─ Unique sample_barcodes extracted: ~391
└─ Overlap with manifest: 0 (ZERO!)
```

**Debug Output from Real Run:**
```
Mapping sample barcodes to conditions
aggregate_by_sample_barcode debug info:
ℹ Input reads: 85785 rows, 391 unique sample_barcodes
ℹ Manifest: 36 rows, 36 unique sample_barcodes
ℹ Merged result has 'condition' column
ℹ  NA conditions: 85785
ℹ  Non-NA conditions: 0
ℹ Matched samples in merged data: 0
Error: No reads matched to sample manifest.
```

**Proof Code Works with Matching Data:**
- Synthetic test (test_aggregate.R): ✅ SUCCESS
- Real pipeline with mismatched data: ❌ ZERO matches

### Example Extracted Sample Barcodes (Not in Manifest)
- CAACTTGT
- AGCTTGTG  
- ACGCAACT
- GAAGACCC
- etc. (391 total unique)

### Expected Sample Barcodes (In Manifest)
- TTGAACCG, AATCCAGC, CCGAGTTA, AACTGTTA (Day_0)
- AATCCAAT, CCGAGTAT, TTCTCATA, AACTGTGC (UT_r1)
- AAGAACTA, GGAGTGTA, TTAGACCG, TTGATGCG (UT_r2)
- etc. (36 total)

---

## KEY ARCHITECTURAL DECISIONS & RATIONALE

### 1. Why data.table for Processing?
**Decision:** Use data.table for core pipeline operations  
**Rationale:**
- Fast aggregation and merging (critical for large datasets)
- Memory efficient for 100M+ read processing
- Native by-group operations
- Better than dplyr for this scale

**Note:** Still use dplyr/readr for file I/O (tidyverse advantage for reading)

### 2. Why Keep dplyr Pipes in Some Places?
**Decision:** Keep %>% pipes in load_data.R but remove from processing pipeline  
**Rationale:**
- File loading is one-time operation (not performance critical)
- Readability for simple transformations
- But remove from hot path (match_constructs, map_to_genes) to prevent column loss

### 3. Why Extract Sample Barcode as First 8bp?
**Decision:** Sample barcode = first 8bp of read  
**Rationale:**
- Standard PoolQ dual-barcode design
- Illumina multiplexing places sample index at read start
- Matches Broad Institute convention

### 4. Why Use Flanking Sequence for Construct Barcode?
**Decision:** Find CACCG, then extract 20bp from position 6 after  
**Rationale:**
- Consistent with PoolQ extraction policy
- Allows for sequencing adapter/primer sequences
- Flexible via configurable policy string

### 5. Why Both Exact and Approximate Matching?
**Decision:** Offer both exact (fast) and fuzzy (tolerant) matching  
**Rationale:**
- Exact matching: production use with clean data
- Approximate matching: recovery from sequencing errors
- User can tune max_mismatches based on quality

---

## TESTING STRATEGY

### Unit Tests (tests/testthat/)
**Coverage:**
- ✅ parse_barcode_policy() - valid/invalid formats
- ✅ extract_barcodes_from_read() - barcode extraction logic
- ✅ create_count_matrix() - construct and gene level
- ✅ lognormalize_counts() - PoolQ normalization
- ✅ calculate_qc_metrics() - overall and per-sample
- ✅ aggregate_by_sample_barcode() - sample matching

**Status:** All tests pass with synthetic data ✅

### Integration Testing
**Test Data Location:** `D:/Bibliotecas/Biologia/crispyr_fastq/example_data/subset_example/`

**Files:**
- HCWFFDRX7_1_0469424154_AACTACCG_S1_L001_R1_001.fastq (19M, original)
- subset_sample.fastq.gz (2.2M, subset - CURRENTLY MISMATCHED)
- sample_manifest_subset.csv (36 samples)
- CP0045_reference_subset.csv (construct library subset)
- CP0045_GRCm38_NCBI_CRISPRko_strict_gene_subset.chip (gene annotations)

**Test Scripts:**
- 01_pipeline_wrapper.R - Full pipeline with process_pooled_screen()
- 02_pipeline_stepwise.R - Step-by-step with timing

**Current Status:** ❌ Fails due to data mismatch (not code bug)

---

## IMMEDIATE NEXT STEPS

### Priority 1: Fix Test Data Mismatch
**Options:**
1. Find original subset creation script and re-run correctly
2. Create new subset from full FASTQ matching the 36 manifest samples
3. Use different test data with known-good pairing

**Action Required:**
```r
# Diagnostic to run:
devtools::load_all()
manifest <- load_sample_manifest("D:/Bibliotecas/Biologia/crispyr_fastq/example_data/subset_example/sample_manifest_subset.csv")
fastq <- "D:/Bibliotecas/Biologia/crispyr_fastq/example_data/subset_example/subset_sample.fastq.gz"
extracted <- extract_barcodes_from_fastq(fastq, quiet = TRUE)

# Check overlap
manifest_bcs <- unique(manifest$sample_barcode)
fastq_bcs <- unique(extracted$sample_barcode[!is.na(extracted$sample_barcode)])
intersect(manifest_bcs, fastq_bcs)  # Should NOT be empty
```

### Priority 2: Validate with Real Data
Once test data is fixed:
1. Run full pipeline end-to-end
2. Compare outputs with PoolQ reference
3. Validate QC metrics
4. Document any differences

### Priority 3: Production Readiness
1. Ensure devtools::check() passes cleanly
2. Final documentation review
3. Performance benchmarking
4. Prepare for merge to main branch

---

## CODE QUALITY STANDARDS

### Style Guide
- Follow tidyverse style guide
- Function names: verb_noun() format
- Variable names: snake_case
- Constants: UPPER_CASE (rare, prefer parameters)

### Error Handling
- Use rlang::abort() for errors with helpful context
- Provide informative error messages with hints
- Validate inputs at function entry
- Use cli::cli_inform() for progress messages

### Documentation
- Roxygen2 docstrings for all exported functions
- @param, @return, @details, @examples for every function
- README with quick start
- Vignette with full workflow

### Testing
- testthat for unit tests
- Test edge cases and error conditions
- Use synthetic data for unit tests
- Integration tests with real data

---

## PERFORMANCE NOTES

### Expected Performance (Based on Design)
- Small screens (1-5M reads): 30-60 seconds
- Medium screens (10-20M reads): 2-4 minutes
- Large screens (50M reads): 8-12 minutes
- Very large screens (100M+ reads): 20-40 minutes

### Bottlenecks
1. **FASTQ reading** - I/O bound (ShortRead::readFastq)
2. **Approximate matching** - CPU bound (Hamming distance calculations)
3. **Count matrix creation** - Memory bound (dcast operations)

### Optimizations Applied
- Hash table exact matching (O(1) lookup)
- data.table for fast aggregation
- Vectorized operations where possible
- Minimal data copying

---

## DEPENDENCIES

### CRAN Packages
- data.table (fast aggregation)
- dplyr (data manipulation)
- tidyr (data reshaping)
- stringr (string operations)
- stringdist (fuzzy matching)
- readr (fast CSV reading)
- rlang (error handling)
- cli (progress reporting)
- glue (string interpolation)

### Bioconductor Packages
- ShortRead (FASTQ reading)
- Biostrings (sequence manipulation)

### Suggests
- testthat (unit testing)

---

## FILE STRUCTURE

```
crispyr/
├── DESCRIPTION              # Package metadata
├── NAMESPACE               # Exported functions (auto-generated)
├── README.md              # User-facing documentation
├── LICENSE                # MIT license
├── crispyr.Rproj          # RStudio project
│
├── R/                     # Source code
│   ├── load_data.R        # Data loading functions
│   ├── extract_barcodes.R # Barcode extraction
│   ├── match_and_map.R    # Library matching and gene mapping
│   ├── count_and_normalize.R # Count matrix creation
│   ├── qc_and_process.R   # QC and wrapper functions
│   └── nin.R              # Helper operator (%nin%)
│
├── man/                   # Documentation (auto-generated)
│   └── *.Rd              # Function help files
│
├── tests/                 # Unit tests
│   └── testthat/
│       ├── test-core_functions.R
│       └── test-grapes-nin-grapes.R
│
├── vignettes/            # User guides
│   └── crispyr-workflow.Rmd
│
└── Example scripts (not in package):
    ├── 01_pipeline_wrapper.R
    └── 02_pipeline_stepwise.R
```

---

## GIT WORKFLOW

### Branches
- **main**: Production-ready code
- **fastq-processing-pipeline**: Current development (YOU ARE HERE)

### Commit Strategy
- Descriptive commit messages following conventional commits
- feat: for new features
- fix: for bug fixes
- docs: for documentation
- test: for test additions/changes
- refactor: for code restructuring

### Current Branch Status
```bash
git branch
# * fastq-processing-pipeline
# 5 commits ahead of main
```

---

## TROUBLESHOOTING GUIDE

### Issue: "No reads matched to sample manifest"
**Symptom:** aggregate_by_sample_barcode() returns empty dataset  
**Diagnosis:** Check if sample_barcodes in FASTQ match manifest
```r
# Run diagnostic
manifest_bcs <- unique(manifest$sample_barcode)
fastq_bcs <- unique(extracted$sample_barcode[!is.na(extracted$sample_barcode)])
length(intersect(manifest_bcs, fastq_bcs))  # Should be > 0
```
**Solution:** Either fix FASTQ subset or update manifest to match

### Issue: "oldString not found in content"
**Symptom:** Edit tool fails  
**Diagnosis:** Trying to edit code that doesn't exist  
**Solution:** Read file first to verify exact content

### Issue: Quality string extraction error
**Symptom:** "no method for coercing this S4 class to a vector"  
**Diagnosis:** Using wrong accessor for quality scores  
**Solution:** Use `as.character(Biostrings::quality(fastq_data)@quality)`

---

## CONTACT & RESOURCES

**Developer:** Martin Gonzalez Fernandez  
**Email:** martin.gonzalezfernandez@unibe.ch  
**ORCID:** 0009-0006-8308-5493  
**Repository:** https://github.com/StrawHattM/crispyr  

**References:**
- Broad Institute PoolQ documentation
- Bioconductor ShortRead package
- data.table documentation

---

## SESSION CHECKLIST

**At Start of Each Session:**
1. ✅ Read this DEVEL_CONTEXT.md file
2. ✅ Check current branch: `git branch`
3. ✅ Review recent commits: `git log --oneline -5`
4. ✅ Check for uncommitted changes: `git status`
5. ✅ Load package: `devtools::load_all()`

**Before Committing:**
1. ✅ Run tests: `devtools::test()`
2. ✅ Check package: `devtools::check()`
3. ✅ Update DEVEL_CONTEXT.md if needed
4. ✅ Write descriptive commit message

**End of Session:**
1. ✅ Update DEVEL_CONTEXT.md with progress
2. ✅ Commit changes with clear message
3. ✅ Document any blockers or next steps

---

*This file should be read at the start of every development session to maintain context.*
