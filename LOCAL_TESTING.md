# Local Testing Guide for crispyr Package

## Prerequisites

Make sure you have installed in R:
```r
install.packages(c("devtools", "roxygen2", "testthat"))
```

## Steps to Run devtools::check()

### 1. Open RStudio
- Navigate to: `D:/Bibliotecas/Biologia/crispyr_fastq/crispyr_repo/crispyr.Rproj`
- This will open the project in RStudio

### 2. Activate the fastq-processing-pipeline branch
```r
# In RStudio console or Terminal tab:
system('git checkout fastq-processing-pipeline')
```

### 3. Generate Documentation (REQUIRED)
This must be done first before running check():
```r
roxygen2::roxygenise()
```

This will create/update:
- `NAMESPACE` - Auto-generated function exports
- `man/*.Rd` - All function documentation files

### 4. Run Package Check
```r
devtools::check()
```

### 5. Expected Output Structure
The check should produce:
- ✓ Loading source code
- ✓ Documenting with roxygen2
- ✓ Checking R code
- ✓ Running tests
- ✓ Building vignettes
- ✓ Checking examples
- ✓ And more...

## What devtools::check() Does

1. **Loads all functions** from R/ directory
2. **Checks NAMESPACE** exports
3. **Validates package metadata** (DESCRIPTION)
4. **Runs all unit tests** (tests/testthat/)
5. **Builds vignettes** (vignettes/)
6. **Checks code style** and common issues
7. **Validates documentation** completeness
8. **Tests function examples**

## Addressing Common Issues

### Issue: "object 'ShortRead' not found"
**Solution**: Install Bioconductor packages
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
```

### Issue: "Missing dependencies"
**Solution**: Install from DESCRIPTION
```r
devtools::install_deps()
```

### Issue: "Test failures"
**Solution**: Run tests individually
```r
devtools::test()
# Or specific test file:
devtools::test_file("tests/testthat/test-core_functions.R")
```

### Issue: "Warning: NAMESPACE not found"
**Solution**: Generate it first
```r
roxygen2::roxygenise()
```

## Quick Checklist

- [ ] Switched to `fastq-processing-pipeline` branch
- [ ] Installed Bioconductor packages (ShortRead, Biostrings)
- [ ] Ran `roxygen2::roxygenise()`
- [ ] Ran `devtools::check()`
- [ ] All tests pass
- [ ] No major warnings/errors

## Loading the Package for Testing

After roxygen2::roxygenise(), you can load and test functions:

```r
# Load all functions
devtools::load_all()

# Test a function
ref_lib <- load_reference_library(
  "example_data/CP0045_Brie_AiO/CP0045_reference_20160120.csv"
)
head(ref_lib)
```

## Running Individual Tests

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-core_functions.R")

# Run with verbose output
devtools::test(reporter = "progress")
```

## Expected Check Output

A successful check should show:
```
-- R CMD check results --------
0 errors ✓ | 0 warnings ✓ | 0 notes ✓
```

Minor notes about "Possibly misspelled words" in documentation are typically acceptable.

## Troubleshooting Resources

- Package development guide: https://r-pkgs.org/
- devtools documentation: https://devtools.r-lib.org/
- testthat guide: https://testthat.r-lib.org/

## Next After Check Passes

Once `devtools::check()` passes cleanly:

1. Install locally: `devtools::install()`
2. Use package: `library(crispyr)`
3. Test with real data (CP0045 example)
4. Compare outputs with PoolQ results
5. Consider pushing to GitHub

---

**Note**: The package is fully implemented. The check will validate that everything is correctly structured and documented for CRAN/Bioconductor submission.
