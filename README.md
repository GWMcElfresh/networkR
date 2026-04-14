# networkR

`networkR` is an R package scaffold for configurable Bayesian network analysis in
factorial and repeated-measures studies. 

I intend to extend this to gaussian copulas/other graphical models, but for now, just Bayesian networks. 


## Current scope

- YAML or programmatic configuration
- `wide_by_tissue` and long-format `panel` data ingestion, with automatic measurement column inference
- Derived design variables
- MICE predictive mean matching imputation
- Hartemink or factor-based discretization with fallback handling
- Tier-aware blacklist construction
- Bootstrap Bayesian network fitting — `stagedtrees` backend (default) and `bnlearn` backend
- `ggraph(layout = "fr")` network plotting; event tree / CEG visualization for the `stagedtrees` backend
- MCP server for natural-language-driven configuration and pipeline execution

## Quick start — MCP-driven workflow (recommended)

The package ships with an MCP server at `inst/mcp_server.R`.  Open the repository
in VS Code; the server starts automatically from `.vscode/mcp.json`.

```
# In GitHub Copilot Chat (or any MCP-compatible client)

1. "Explore my data directory at data/"
   → networkr_explore_data_directory reports column roles and a reviewScaffold per file.

2. "Generate a config from data/subjectIdTable_TB.csv"
   → networkr_generate_configuration_from_data returns a YAML draft and a columnRoles scaffold.

3. Review the scaffold. Move columns between roles as needed, then:
   "Refine the config using this corrected scaffold: <paste edited columnRoles block>"
   → networkr_refine_configuration_from_scaffold regenerates YAML from the corrected assignments.

4. "Save the config to configs/my_analysis.yml"
   → networkr_write_configuration_file writes and immediately validates the file.

5. "Run the pipeline from configs/my_analysis.yml"
   → networkr_run_pipeline_from_config runs the full analysis.

6. "Describe the results"
   → networkr_describe_last_pipeline_result summarizes the most recent run.
```

See `vignettes/networkR-mcp.Rmd` for a full step-by-step walkthrough.

## Quick start — programmatic

```r
library(networkR)

configuration <- ParseConfiguration("configs/my_analysis.yml")
analysisResult <- RunPipeline(configuration)

analysisResult$plotObjects$fullModel
```

`configs/lung_ldln_bidirectional_template.yml` is an annotated reference config for
`wide_by_tissue` bidirectional analyses; copy and edit the column-name fields to match
your data.

## Status

This repository contains the package foundation. Remaining work includes joint
paired-tissue modeling, mediation workflows, and correlation reporting.
