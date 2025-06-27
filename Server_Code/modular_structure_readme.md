# Modular Shiny Server Structure

This directory contains the modularized server code for the miRQuest Shiny application. The code has been separated into logical, functional modules while maintaining all reactive variables and their relationships.

## File Structure

### Core Files
- `server.R` - Main server file with reactive variable definitions and file sourcing
- `global.R` - Libraries, settings, and helper function sourcing

### Backend Server Code Modules (`Backend_Server_Code/`)

1. **`data_input_handlers.R`**
   - Module 0A: miRNA counts data reading and processing
   - Module 0B: Metadata file handling and UI updates

2. **`data_subsetting.R`**
   - Module 2: Data filtering and subsetting functionality
   - Dynamic UI generation for filter values
   - Download handlers for subsetted data

3. **`exploratory_visualization.R`**
   - Module 1: Exploratory data visualization
   - Stacked column charts, MDS plots, heatmaps
   - Data normalization and transformation

4. **`differential_mirna_analysis.R`**
   - Module 3: DESeq2 differential miRNA expression analysis
   - Volcano plots, significance testing
   - Alternative limma::voom analysis
   - Boxplot generation for individual miRNAs

5. **`differential_gene_analysis.R`**
   - Module 6: DESeq2 differential gene expression analysis
   - mRNA data processing and analysis
   - Results filtering and visualization

6. **`correlation_analysis.R`**
   - Module 7: miRNA-mRNA correlation analysis
   - Negative correlation identification
   - Database support validation
   - Correlation plotting and filtering

7. **`target_prediction_all.R`**
   - Module 4: Bulk target gene prediction for all DE miRNAs
   - Pathway overrepresentation analysis
   - Up/down-regulated miRNA separation
   - Progress tracking for large analyses

8. **`target_prediction_single.R`**
   - Module 5: Single miRNA target prediction
   - Individual pathway analysis
   - Dotplot and barplot generation

9. **`pathway_visualization.R`**
   - Chord plot generation
   - Network plot visualization
   - Interactive pathway-miRNA relationships

10. **`future_expansion_code.R`**
    - Commented-out code for future development
    - Alternative implementations
    - Experimental features

## Reactive Variables

All reactive variables are defined in the main `server.R` file and are accessible across all modules:

### Species and Database Reactives
- `mart()` - Biomart connection
- `org_string()` - Organism-specific gene symbols
- `org_abbrev()` - Species abbreviation
- `org_db()` - Organism database
- `organism_reactive()` - Full organism name

### Data Storage Reactives
- `longData()`, `wideData()` - miRNA count data
- `metadata()` - Sample metadata
- `DESEQ_obj()` - DESeq2 object
- `res_significant()` - Significant DE miRNAs
- `mrna_res_significant()` - Significant DE genes

### Analysis Results Reactives
- `predicted_target_reactive()` - All predicted targets
- `neg_cor_reactive()` - Negative correlations
- `All_miRNA_Pathways_reactive()` - Pathway analysis results
- And many more...

## Key Features

### Tidyverse Code Style
- Proper indentation and spacing
- Function arguments on separate lines for clarity
- Consistent naming conventions

### Error Handling
- Notification systems for long-running processes
- Validation checks for required inputs
- Graceful handling of missing data

### Download Capabilities
- CSV downloads for all major results
- Plot downloads in multiple formats
- Consistent file naming with timestamps

### Progress Tracking
- Progress bars for intensive computations
- User feedback during database queries
- Notification management

## Usage

The modular structure is automatically sourced in `server.R`:

```r
files_to_source <- list.files("Backend_Server_Code")

for (file in files_to_source) {
    source(paste0("Backend_Server_Code/", file), local = TRUE)
}
```

This approach maintains all reactive relationships while providing clean separation of concerns. Each module focuses on a specific functionality while having access to all shared reactive variables.

## Benefits

1. **Maintainability** - Each module handles specific functionality
2. **Readability** - Code is organized by purpose
3. **Debugging** - Issues can be isolated to specific modules
4. **Collaboration** - Multiple developers can work on different modules
5. **Testing** - Individual modules can be tested independently
6. **Scalability** - New features can be added as new modules

## Dependencies

Each module assumes access to:
- All reactive variables defined in `server.R`
- All libraries loaded in `global.R`
- Helper functions sourced from the `scripts/` directory
- Standard Shiny server environment (`input`, `output`, `session`)
