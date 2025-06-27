# Create_Exe_File.R
# This script creates a Windows executable file for launching the miRQuest Shiny app

# Check if running on Windows
if (.Platform$OS.type != "windows") {
  stop("This script is only compatible with Windows. For non-Windows users, please use the .Rproj file method.")
}

# Check if shiny.exe is installed, if not, install it
if (!require("shiny.exe", quietly = TRUE)) {
  message("Installing shiny.exe package...")
  install.packages("shiny.exe")
  library(shiny.exe)
}

# Create the executable file
shiny.exe()

message("miRQuest.bat file has been created successfully!")
message("You can now double-click the miRQuest.bat file to launch the application.")