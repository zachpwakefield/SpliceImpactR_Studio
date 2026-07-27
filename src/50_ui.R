sir_section_card <- function(title, ..., class = NULL, controls = NULL) {
  div(
    class = paste("panel-card", class %||% ""),
    div(
      class = "panel-card__header",
      div(class = "panel-card__title", title),
      if (!is.null(controls)) div(class = "panel-card__controls", controls)
    ),
    div(class = "panel-card__body", ...)
  )
}

sir_download_row <- function(...) {
  div(class = "download-row", ...)
}

sir_help_text <- function(...) {
  div(class = "sir-help help-note", ...)
}

sir_table_block <- function(title, output_id, caption_id, search_id = NULL) {
  tagList(
    div(
      class = "sir-table-toolbar",
      div(class = "sir-table-title", title),
      if (!is.null(search_id)) {
        textInput(
          search_id,
          label = NULL,
          placeholder = "Filter this preview"
        )
      }
    ),
    uiOutput(caption_id, class = "sir-table-caption"),
    DT::DTOutput(output_id)
  )
}

sir_action_button <- function(id, label, primary = FALSE, disabled = FALSE) {
  actionButton(
    id,
    label,
    class = paste(
      "btn",
      if (isTRUE(primary)) "btn-primary" else "btn-default"
    ),
    disabled = if (isTRUE(disabled)) "disabled" else NULL
  )
}

sir_ui <- function() {
  fluidPage(
    theme = bslib::bs_theme(
      version = 5,
      bg = "#f4efe6",
      fg = "#1d2926",
      primary = "#16645a",
      secondary = "#b25c33",
      success = "#3c7a55",
      info = "#4c8e9a",
      warning = "#bd7a2d",
      danger = "#a94b4b",
      base_font = bslib::font_collection(
        "Avenir Next",
        "Segoe UI",
        "sans-serif"
      ),
      heading_font = bslib::font_collection(
        "Iowan Old Style",
        "Georgia",
        "serif"
      )
    ),
    tags$head(
      tags$title("SpliceImpactR Studio"),
      tags$meta(
        name = "description",
        content = "Public analysis workspace for SpliceImpactR"
      ),
      includeCSS(file.path(SIR_APP_DIR, "www", "spliceimpactr.css")),
      tags$script(HTML(
        "
        Shiny.addCustomMessageHandler('sir-button-state', function(message) {
          var element = document.getElementById(message.id);
          if (!element) return;
          element.disabled = Boolean(message.disabled);
          element.setAttribute('aria-disabled', message.disabled ? 'true' : 'false');
          element.classList.toggle('sir-disabled', Boolean(message.disabled));
        });
        "
      ))
    ),
    div(
      class = "app-shell",
      div(
        class = "masthead",
        div(class = "masthead__eyebrow", "Shiny Workspace For SpliceImpactR · v4"),
        h1("SpliceImpactR Studio"),
        p(
          "Load annotations, sample manifests, user-defined event tables, protein features, ",
          "and PPI data; then run differential inclusion, transcript mapping, sequence/frame ",
          "comparison, domain and interaction analysis, transcript browsing, and gene-set enrichment ",
          "from one coherent interface."
        ),
        div(
          class = "masthead__chips",
          div(class = "masthead__chip", "Sample-based rMATS + HIT Index workflow"),
          div(class = "masthead__chip", "User-supplied raw and post-DI imports"),
          div(class = "masthead__chip", "Sequence, domain, PPI, and enrichment layers"),
          div(
            class = "masthead__chip",
            "Bundled full reference + browser-safe uploads"
          )
        ),
        uiOutput("runtime_banner")
      ),
      fluidRow(
        column(
          width = 4,
          div(
            class = "control-rail",
            bslib::accordion(
              id = "control_sections",
              open = c("Quick Start", "References"),
              bslib::accordion_panel(
                "Quick Start",
                p(
                  class = "help-note",
                  "Demo mode loads example annotations, a bundled manifest, ",
                  "reference protein features, and sample-level events. ",
                  "Run Full Pipeline loads the selected PPI network on demand."
                ),
                sir_action_button(
                  "load_demo_workspace",
                  "Load Demo Workspace",
                  primary = TRUE
                ),
                sir_action_button(
                  "run_all_pipeline",
                  "Run Full Pipeline",
                  disabled = TRUE
                ),
                sir_action_button("clear_analysis", "Clear Results"),
                sir_action_button("new_workspace", "New Workspace")
              ),
              bslib::accordion_panel(
                "References",
                radioButtons(
                  "annotation_mode",
                  "Annotation source",
                  choices = sir_reference_mode_choices(),
                  selected = sir_default_reference_mode()
                ),
                uiOutput("reference_source_controls"),
                uiOutput("reference_description"),
                uiOutput("reference_storage_status"),
                sir_action_button(
                  "load_reference",
                  "Load Annotations",
                  primary = TRUE
                )
              ),
              bslib::accordion_panel(
                "Event Inputs",
                radioButtons(
                  "input_mode",
                  "Evidence source",
                  choices = c(
                    "Bundled sample manifest" = "demo",
                    "Self-contained rMATS project ZIP" = "project",
                    "Normalized raw-event table" = "raw",
                    "Post-DI event table" = "di"
                  ),
                  selected = "demo"
                ),
                uiOutput("input_mode_controls"),
                sir_action_button("load_input", "Load Event Data", primary = TRUE),
                tags$hr(),
                uiOutput("condition_controls")
              ),
              bslib::accordion_panel(
                "Protein Features",
                radioButtons(
                  "feature_mode",
                  "Feature source",
                  choices = c(
                    "Packaged test subset" = "demo",
                    "Query Ensembl / BioMart" = "query",
                    "Upload feature table" = "upload"
                  ),
                  selected = "demo"
                ),
                conditionalPanel(
                  "input.feature_mode === 'query'",
                  numericInput(
                    "biomart_release",
                    "BioMart/Ensembl release",
                    value = 111,
                    min = 1,
                    step = 1
                  ),
                  checkboxGroupInput(
                    "feature_sources",
                    "Feature databases",
                    choices = SIR_FEATURE_SOURCES,
                    selected = c("interpro", "pfam", "signalp")
                  ),
                  checkboxInput(
                    "combine_feature_overlaps",
                    "Combine overlapping instances",
                    value = FALSE
                  ),
                  p(
                    class = "help-note",
                    "Queries are cached only in the configured managed feature cache."
                  )
                ),
                conditionalPanel(
                  "input.feature_mode === 'upload'",
                  fileInput(
                    "feature_upload",
                    "Feature table",
                    accept = c(
                      ".csv", ".tsv", ".txt",
                      ".csv.gz", ".tsv.gz", ".txt.gz"
                    )
                  ),
                  radioButtons(
                    "feature_upload_mode",
                    "Uploaded schema",
                    choices = c(
                      "Normalized SpliceImpactR features" = "normalized",
                      "Manual protein coordinates" = "manual"
                    ),
                    selected = "normalized"
                  )
                ),
                sir_action_button(
                  "load_features",
                  "Load Protein Features",
                  primary = TRUE,
                  disabled = TRUE
                ),
                uiOutput("feature_coverage")
              ),
              bslib::accordion_panel(
                "PPI + Background",
                radioButtons(
                  "ppi_mode",
                  "PPI source",
                  choices = c(
                    "Bundled SpliceImpactR network" = "default",
                    "Upload normalized network" = "upload"
                  ),
                  selected = "default"
                ),
                conditionalPanel(
                  "input.ppi_mode === 'upload'",
                  fileInput(
                    "ppi_upload",
                    "PPI table",
                    accept = c(
                      ".csv", ".tsv", ".txt",
                      ".csv.gz", ".tsv.gz", ".txt.gz"
                    )
                  )
                ),
                sir_action_button("load_ppi", "Load PPI Interactions"),
                uiOutput("ppi_status"),
                tags$hr(),
                radioButtons(
                  "background_source",
                  "Enrichment background",
                  choices = c(
                    "All annotated feature-bearing genes" = "annotated",
                    "Genes observable in the rMATS project" = "hit_index",
                    "My transcript IDs" = "user-given"
                  ),
                  selected = "annotated"
                ),
                conditionalPanel(
                  "input.background_source === 'user-given'",
                  textAreaInput(
                    "background_transcripts",
                    "Transcript IDs",
                    rows = 5,
                    placeholder = "One ENST ID per line or comma-separated"
                  )
                ),
                sir_action_button(
                  "build_background",
                  "Build Enrichment Background",
                  disabled = TRUE
                ),
                uiOutput("background_status")
              ),
              bslib::accordion_panel(
                "Analysis Controls",
                numericInput(
                  "min_reads",
                  "Minimum total reads",
                  value = 10,
                  min = 0,
                  step = 1
                ),
                sliderInput(
                  "min_prop",
                  "Minimum sample proportion",
                  min = 0,
                  max = 1,
                  value = 0.5,
                  step = 0.05
                ),
                p(
                  class = "help-note",
                  "Missing AFE/ALE sites use the maximum observed depth within the same event."
                ),
                textInput("cooks_cutoff", "Cook's distance cutoff", value = "Inf"),
                numericInput(
                  "padj_threshold",
                  "FDR cutoff",
                  value = 0.05,
                  min = 0,
                  max = 1,
                  step = 0.01
                ),
                numericInput(
                  "dpsi_threshold",
                  "Absolute delta PSI cutoff",
                  value = 0.10,
                  min = 0,
                  max = 1,
                  step = 0.01
                ),
                selectInput(
                  "pair_mode",
                  "Transcript pairing mode",
                  choices = c(
                    "All compatible isoforms" = "multi",
                    "Exactly one INC and one EXC row" = "paired"
                  ),
                  selected = "multi"
                ),
                checkboxInput(
                  "show_protein_domains",
                  "Keep whole-protein domain detail",
                  value = FALSE
                ),
                sir_action_button(
                  "run_di",
                  "Run Differential Inclusion",
                  disabled = TRUE
                ),
                sir_action_button(
                  "map_events",
                  "Map Significant Events",
                  disabled = TRUE
                ),
                sir_action_button(
                  "run_consequences",
                  "Run Sequence / Domain / PPI",
                  primary = TRUE,
                  disabled = TRUE
                )
              ),
              bslib::accordion_panel(
                "Filters",
                textInput(
                  "input_search",
                  "Workspace input filter",
                  placeholder = "Gene, event, sample, or condition"
                ),
                textInput(
                  "analysis_search",
                  "Analysis result filter",
                  placeholder = "Gene, transcript, protein, or event"
                ),
                selectInput(
                  "preview_rows",
                  "Rows per preview",
                  choices = c(10, 25, 50, 100, 250),
                  selected = 25
                )
              )
            )
          )
        ),
        column(
          width = 8,
          uiOutput("metric_cards"),
          div(
            class = "panel-shell",
            tabsetPanel(
              id = "main_tabs",
              tabPanel(
                "Workspace",
                sir_section_card("Workspace status", uiOutput("workspace_status")),
                tabsetPanel(
                  id = "workspace_tabs",
                  tabPanel(
                    "Input",
                    sir_section_card(
                      "Current sample manifest or event input",
                      sir_table_block(
                        "Input preview",
                        "input_preview",
                        "input_preview_caption"
                      )
                    )
                  ),
                  tabPanel(
                    "Annotations",
                    sir_section_card(
                      "Annotation summary and preview",
                      sir_table_block(
                        "Current annotation rows",
                        "reference_table",
                        "reference_table_caption",
                        "reference_search"
                      )
                    )
                  ),
                  tabPanel(
                    "Protein features",
                    sir_section_card(
                      "Protein feature preview",
                      sir_table_block(
                        "Current feature rows",
                        "feature_table",
                        "feature_table_caption",
                        "feature_search"
                      )
                    )
                  ),
                  tabPanel(
                    "PPI / background",
                    sir_section_card(
                      "PPI and enrichment background",
                      sir_table_block(
                        "Current resource rows",
                        "resource_table",
                        "resource_table_caption",
                        "resource_search"
                      )
                    )
                  )
                )
              ),
              tabPanel(
                "DI",
                sir_section_card(
                  "Differential inclusion overview",
                  plotOutput("di_plot", height = 360),
                  uiOutput("di_plot_note")
                ),
                sir_section_card(
                  "Significant events",
                  sir_table_block(
                    "Current significant event-form rows",
                    "di_table",
                    "di_table_caption"
                  )
                )
              ),
              tabPanel(
                "Consequences",
                uiOutput("analysis_freshness"),
                uiOutput("analysis_message"),
                uiOutput("consequence_message"),
                tabsetPanel(
                  id = "consequence_tabs",
                  tabPanel(
                    "Mapping",
                    sir_section_card("Matched event summary", uiOutput("mapping_summary")),
                    sir_section_card(
                      "Matched event-to-transcript table",
                      sir_table_block(
                        "Matched rows",
                        "matched_table",
                        "matched_table_caption"
                      )
                    )
                  ),
                  tabPanel(
                    "Sequence",
                    sir_section_card(
                      "Alignment summary",
                      plotOutput("alignment_plot", height = 340)
                    ),
                    sir_section_card(
                      "Length comparison",
                      plotOutput("length_plot", height = 340)
                    ),
                    uiOutput("sequence_note"),
                    sir_section_card(
                      "Sequence/frame table",
                      sir_table_block(
                        "Transcript-pair sequence results",
                        "sequence_table",
                        "sequence_table_caption"
                      )
                    )
                  ),
                  tabPanel(
                    "Domains",
                    sir_section_card(
                      "Domain change summary",
                      plotOutput("domain_plot", height = 360),
                      uiOutput("domain_note")
                    ),
                    sir_section_card(
                      "Domain-level event table",
                      sir_table_block(
                        "Domain changes",
                        "domain_table",
                        "domain_table_caption",
                        "domain_search"
                      )
                    )
                  ),
                  tabPanel(
                    "PPI",
                    sir_section_card(
                      "PPI summary",
                      plotOutput("ppi_plot", height = 360),
                      uiOutput("ppi_note")
                    ),
                    sir_section_card(
                      "PPI event table",
                      sir_table_block(
                        "PPI changes",
                        "ppi_table",
                        "ppi_table_caption",
                        "ppi_search"
                      )
                    )
                  ),
                  tabPanel(
                    "Integrated",
                    sir_section_card(
                      "Integrated event summary",
                      plotOutput("integrated_plot", height = 620),
                      uiOutput("integrated_summary")
                    )
                  )
                )
              ),
              tabPanel(
                "Explore",
                tabsetPanel(
                  id = "explore_tabs",
                  tabPanel(
                    "Event probe",
                    sir_section_card(
                      "Probe one event across samples",
                      fluidRow(
                        column(
                          4,
                          selectizeInput(
                            "probe_gene",
                            "Gene",
                            choices = character(),
                            options = list(placeholder = "Select gene")
                          ),
                          uiOutput("probe_event_type_ui"),
                          uiOutput("probe_event_ui"),
                          sir_action_button(
                            "run_probe",
                            "Probe Event",
                            primary = TRUE,
                            disabled = TRUE
                          )
                        ),
                        column(
                          8,
                          plotOutput("probe_plot", height = 320),
                          uiOutput("probe_note")
                        )
                      )
                    )
                  ),
                  tabPanel(
                    "Transcript browser",
                    sir_section_card(
                      "Direct transcript-to-transcript comparison",
                      fluidRow(
                        column(
                          4,
                          selectizeInput(
                            "browser_gene",
                            "Gene",
                            choices = character(),
                            options = list(placeholder = "Select gene")
                          ),
                          selectInput(
                            "transcript_a",
                            "Reference transcript",
                            choices = character()
                          ),
                          selectInput(
                            "transcript_b",
                            "Comparison transcript",
                            choices = character()
                          ),
                          radioButtons(
                            "browser_view",
                            "View",
                            choices = c(
                              "Transcript" = "transcript",
                              "Protein" = "protein"
                            ),
                            selected = "transcript",
                            inline = TRUE
                          ),
                          selectizeInput(
                            "browser_feature_databases",
                            "Feature databases",
                            choices = character(),
                            multiple = TRUE
                          ),
                          checkboxInput(
                            "browser_combine_domains",
                            "Combine repeated domains",
                            value = FALSE
                          ),
                          sir_action_button(
                            "run_browser",
                            "Render Transcript Browser",
                            primary = TRUE,
                            disabled = TRUE
                          )
                        ),
                        column(
                          8,
                          plotOutput("browser_plot", height = 420),
                          uiOutput("browser_note")
                        )
                      )
                    )
                  ),
                  tabPanel(
                    "HIT / PSI",
                    sir_section_card(
                      "Global HIT index and PSI summaries",
                      fluidRow(
                        column(
                          4,
                          selectInput(
                            "overview_event_type",
                            "PSI overview event type",
                            choices = sir_event_type_choices(SIR_SUPPORTED_EVENTS),
                            selected = "SE"
                          ),
                          sir_action_button(
                            "run_overview",
                            "Run HIT / PSI Summaries",
                            primary = TRUE,
                            disabled = TRUE
                          ),
                          uiOutput("overview_event_types_note")
                        ),
                        column(
                          8,
                          plotOutput("hit_plot", height = 360),
                          plotOutput("psi_plot", height = 380)
                        )
                      )
                    )
                  )
                )
              ),
              tabPanel(
                "Enrichment",
                sir_section_card(
                  "Gene-set enrichment",
                  fluidRow(
                    column(
                      4,
                      selectInput(
                        "enrichment_mode",
                        "Foreground mode",
                        choices = c(
                          "DI" = "di",
                          "Domain" = "domain",
                          "PPI" = "ppi"
                        ),
                        selected = "di"
                      ),
                      selectInput(
                        "enrichment_gene_id_type",
                        "Gene ID type",
                        choices = c("Ensembl" = "ensembl", "Symbol" = "symbol"),
                        selected = "ensembl"
                      ),
                      checkboxGroupInput(
                        "enrichment_sources",
                        "Sources",
                        choices = SIR_ENRICHMENT_SOURCES,
                        selected = c("GO:BP", "GO:MF", "GO:CC")
                      ),
                      numericInput(
                        "enrichment_min_size",
                        "Minimum gene-set size",
                        value = 10,
                        min = 1
                      ),
                      numericInput(
                        "enrichment_max_size",
                        "Maximum gene-set size",
                        value = 2000,
                        min = 2
                      ),
                      numericInput(
                        "enrichment_fdr",
                        "Adjusted p-value cutoff",
                        value = 0.05,
                        min = 0,
                        max = 1,
                        step = 0.01
                      ),
                      radioButtons(
                        "enrichment_plot_type",
                        "Plot type",
                        choices = c("Dot" = "dot", "Bar" = "bar"),
                        selected = "dot",
                        inline = TRUE
                      ),
                      sir_action_button(
                        "run_enrichment",
                        "Run Enrichment",
                        primary = TRUE,
                        disabled = TRUE
                      )
                    ),
                    column(
                      8,
                      uiOutput("enrichment_context"),
                      plotOutput("enrichment_plot", height = 380),
                      uiOutput("enrichment_note")
                    )
                  ),
                  sir_table_block(
                    "Significant gene sets",
                    "enrichment_table",
                    "enrichment_table_caption",
                    "enrichment_search"
                  )
                )
              ),
              tabPanel(
                "Exports",
                sir_section_card(
                  "Download analysis products",
                  p(
                    class = "help-note",
                    "Exports use complete current in-memory objects; previews remain bounded."
                  ),
                  sir_download_row(
                    downloadButton("download_raw", "Raw events"),
                    downloadButton("download_di", "DI results"),
                    downloadButton("download_sig", "Significant DI"),
                    downloadButton("download_matched", "Matched events")
                  ),
                  sir_download_row(
                    downloadButton("download_sequence", "Sequence results"),
                    downloadButton("download_domains", "Domain changes"),
                    downloadButton("download_ppi_results", "PPI changes"),
                    downloadButton("download_final", "Integrated hits")
                  ),
                  sir_download_row(
                    downloadButton("download_enrichment", "Enrichment"),
                    downloadButton("download_log", "Activity log"),
                    downloadButton("download_provenance", "Provenance JSON")
                  )
                ),
                sir_section_card(
                  "Workspace manifest",
                  p(
                    class = "help-note",
                    "Reference metadata, input source, counts, parameters, and stage signatures."
                  ),
                  verbatimTextOutput("provenance_preview")
                )
              ),
              tabPanel(
                "Logs",
                sir_section_card(
                  "Execution log",
                  verbatimTextOutput("activity_log")
                )
              )
            )
          )
        )
      )
    )
  )
}
