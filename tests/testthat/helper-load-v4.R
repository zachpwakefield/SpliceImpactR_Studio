test_project_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
options(spliceimpactr.app_dir = test_project_root)

source(file.path(test_project_root, "src", "00_config.R"))
source(file.path(test_project_root, "src", "05_resources.R"))
source(file.path(test_project_root, "src", "10_utils.R"))
source(file.path(test_project_root, "src", "20_io.R"))
source(file.path(test_project_root, "src", "30_state.R"))
source(file.path(test_project_root, "src", "40_pipeline.R"))
